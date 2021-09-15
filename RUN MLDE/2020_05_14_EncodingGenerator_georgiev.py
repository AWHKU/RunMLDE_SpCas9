"""
This module contains a class needed to generate all encoding types.
"""
# Import necessary MLDE modules
# Import necessary MLDE modules
from Support.Encode.MolBioInfo import all_aas, allowed_aas
from Support.Encode.GeorgievParams import georgiev_parameters
#from Support.Encode.TapeModelLocations import tape_model_locations

all_aas = ("A", "C", "D", "E", "F", "G", "H", "I", "K", "L",
           "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
allowed_aas = set(all_aas)



# Import other necessary modules
import os
import re
import warnings
import pickle
import subprocess
import numpy as np
from itertools import product
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from time import strftime

# Write a regex that splits protein amino acid indices into amino acid code and
# integer
_aa_ind_splitter = re.compile("([A-Za-z])([0-9]+)")

# Get the directory of this file
_filedir = os.path.dirname(os.path.abspath(__file__))

#===============================================================================
#============================== Helper Functions ===============================
#===============================================================================
# Write a function that normalizes encodings
def _normalize_encodings(unnormalized_encodings):
    """
    Takes a tensor of embeddings, flattens the internal dimensions, then mean-centers
    and unit-scales. After scaling, the matrix is repacked to its original shape.
    
    Parameters
    ----------
    unnormalized_encodings: Numpy array of shape N x A x D, where N is the number
        of combinations in the design space, A is the number of amino acids in 
        the combination, and D is the dimensionality of the base encoding. This
        array contains the unnormalized MLDE encodings.
        
    Returns
    -------
    normalized_encodings: Numpy array with the same shape as the input, only
        with encodings mean-centered and unit-scaled
    """
    # Raise an error if the input array is not 3D
    assert len(unnormalized_encodings.shape) == 3, "Input array must be 3D"
    
    # Get the length of a flattened array
    flat_shape = np.prod(unnormalized_encodings.shape[1:])

    # Flatten the embeddings along the internal axes
    flat_encodings = np.reshape(unnormalized_encodings,
                                 [len(unnormalized_encodings), flat_shape])

    # Mean-center and unit variance scale the embeddings.
    means = flat_encodings.mean(axis=0)
    stds = flat_encodings.std(axis=0)
    normalized_flat_encodings = (flat_encodings - means)/stds
    
    # Reshape the normalized embeddings back to the original form
    normalized_encodings = np.reshape(normalized_flat_encodings,
                                      unnormalized_encodings.shape)
    
    return normalized_encodings

# Define a class for generating embeddings
class EncodingGenerator():
    """
    The class which contains all information needed to generate encodings.
    
    Parameters
    ----------
    encoding: str: Choice of "learned", "georgiev", and "onehot". This dictates
        how the combinatorial space will be encoded.
    protein_name: str: Nickname for the combinatorial space that will be built.
    fasta_path: str (default = None): Path to the fasta file containing the parent
        protein sequence. This argument is required when using learned embeddings;
        it is ignored when using other encodings.
    target_protein_indices: list of str: Positions in the protein to encode. Must
        take the format 'WTaaPos', where 'WTaa' is the wild-type amino acid at
        position 'Pos' (1-indexed): e.g. V20 means that position 20 in the protein
        given by 'fasta_path' has Valine in the wild type, and encodings should be
        built for it. A list of positions defines the combinatorial space to encode.
        This argument is required when using learned embeddings; otherwise it is
        ignored.
    n_positions_combined: int: The number of amino acids to combine. Ignored 
        when using learned embeddings, required otherwise.
    output: str: Location to save all data. By default, this is the current
        working directory.
    
    Returns
    -------
    None. Outputs are saved to the location given by 'output'.
    """
    # Initialize the embedding
    def __init__(self, encoding, protein_name,
                 fasta_path = None, target_protein_indices = None, 
                 n_positions_combined = None, output = os.getcwd()):

        # Assign all inputs as instance variables
        self._encoding = encoding.lower()
        self._fasta_path = fasta_path
        self._target_protein_indices = target_protein_indices
        self._protein_name = protein_name
        self._output = output
        if encoding in {"georgiev", "onehot"}:
            
            # Build output directories
            self._build_output_dirs()        
            # Build the list of combinations for the position and the dictionaries 
            # linking position index to combo
            self._build_combo_dicts()    
            self.generate_encodings()    
        else:
            raise AssertionError("Unknown encoding")
        
        


    # Write a function that builds output directories
    def _build_output_dirs(self):
        """
        Self-explanatory: Build necessary directories for saving data.
        """
        # Get the start time
        init_time = strftime("%Y%m%d-%H%M%S")
        
        # Build the output directories consistent for all encodings
        self._encoding_output = os.path.join(self.output, init_time, "Encodings")
        os.makedirs(self.encoding_output)
        
        # Build the output directories only used only for generating learned embeddings
        self._fasta_output = os.path.join(self.output, init_time, "Fastas")
        os.makedirs(self.fasta_output)

    # Write a function that analyzes the input fasta file
    # input is a fasta files with all combo sequences of interests made in R
    def _make_combos_from_input_fasta(self):
        """
        Loads the input fasta file and makes sure it passes a number of checks.
        Sets the variable self._wt_seq, which contains the sequence in the input
        fasta.
        """
        # Check to make sure the file exists
        if not os.path.exists(self.fasta_path):
            raise IOError("Cannot locate '{}'".format(self.fasta_path))

        # Load the fasta file
        with open(self.fasta_path, "r") as f:            
            # Open the fasta file and extract all sequences
            fasta_seqs = list(SeqIO.parse(f, "fasta"))

        def phrase_fasta_header(fasta_seq_id):
            aas=fasta_seq_id.split("_")[1]
            aas=aas.split("-")
            #aa_pos=[aa[2:-1] for aa in aas]
            mut_aa=tuple([aa[-1] for aa in aas])
            return mut_aa

        combo_list_from_fasta =[]    
        for i in fasta_seqs:
            combo_list_from_fasta.append(phrase_fasta_header(i.id))
        self._all_combos = combo_list_from_fasta 

            
    # Write a function that produces dictionaries linking combo and index in the
    # output encodings
    def _build_combo_dicts(self):
        """
        Builds dictionaries which link the identity of a combination (e.g. ACTV)
        to that combination's index in an encoding array, and vice versa. Both
        dictionaries are saved to disk.
        """
        # Identify all possible combinations
        # change it to make combo for a list of specific variants from fasta header
        self._make_combos_from_input_fasta()
        self._combi_space = len(self._all_combos)
        all_combos = self._all_combos 
        self._n_positions_combined=len(all_combos[1])
        # Link combo to index and index to combo
        combo_to_index = {"".join(combo): i for i, combo in enumerate(self.all_combos)}
        self._index_to_combo = {i: "".join(combo) for i, combo in enumerate(self.all_combos)}
        
        # Save the constructed dictionaries
        with open(os.path.join(self.encoding_output, f"{self.protein_name}_{self.encoding}_ComboToIndex.pkl"), "wb") as f:
            pickle.dump(combo_to_index, f)
        with open(os.path.join(self.encoding_output, f"{self.protein_name}_{self.encoding}_IndexToCombo.pkl"), "wb") as f:
            pickle.dump(self.index_to_combo, f)
        


    # Write a function that generates onehot encodings
    def _generate_onehot(self):
        """
        Builds a onehot encoding for a given combinatorial space.
        """
        # Make a dictionary that links amino acid to index
        one_hot_dict = {aa: i for i, aa in enumerate(all_aas)}
    
        # Build an array of zeros
        onehot_array = np.zeros([len(self.all_combos), self.n_positions_combined, 20])
        
        onehot_array = np.zeros([len(all_combos), n_positions_combined, 20])
        # Loop over all combos. This should all be vectorized at some point.
        for i, combo in enumerate(self.all_combos):
            
            # Loop over the combo and add ones as appropriate
            for j, character in enumerate(combo):
                
                # Add a 1 to the appropriate position
                onehot_ind = one_hot_dict[character]
                onehot_array[i, j, onehot_ind] = 1
                
        # Return the array
        return onehot_array
        
    # Write a function that generates georgiev encodings
    def _generate_georgiev(self):
        """
        Encodes a given combinatorial space with Georgiev parameters.
        """
        # Now build all encodings for every combination
        unnorm_encodings = np.empty([len(self.all_combos), 
                                       self.n_positions_combined, 19])
        for i, combo in enumerate(self.all_combos):
            unnorm_encodings[i] = [[georgiev_param[character] for georgiev_param
                                    in georgiev_parameters] for character in combo]
            
        return unnorm_encodings
        

    #===========================================================================
    #============================== Public Methods =============================
    #===========================================================================    
    # Write a function that generates encodings
    def generate_encodings(self, n_batches = 1):
        """
        Generates encodings based on the self.encoding instance variable of the
        encoding. Also performs KS-sampling if desired. Note that this class is
        not currently set up to be reused. In other words, a new class should
        be instantiated for generating a new set of encodings. 
        
        Parameters
        ----------
        n_batches: int: The number of batches to split the job into. TAPE heavily
            uses system RAM, and splitting into batches lowers the memory 
            requirements.
            
        Returns
        -------
        None. Will save normalized and unnormalized encodings to disk. If KS-
            sampling is performed, this will be saved to disk as well.
        """
        # Generate the appropriate encoding
        if self.encoding == "georgiev":
            unnormalized_embeddings = self._generate_georgiev()
        elif self.encoding == "onehot":
            
            # Get the embeddings
            onehot_array = self._generate_onehot()
            
            # Save the encodings
            savename = os.path.join(self.encoding_output,
                                    f"{self.protein_name}_onehot_UnNormalized.npy")
            np.save(savename, onehot_array)
            
            # Return
            return None
            
        else:
            raise ValueError(f"Unknown encoding type {self.encoding}")
        
        # Normalize embeddings
        # Reshape the normalized embeddings back to the original form
        normalized_embeddings = _normalize_encodings(unnormalized_embeddings)

        # Create filenames for saving the embeddings
        unnormalized_savename = os.path.join(self.encoding_output, 
                                             f"{self.protein_name}_{self.encoding}_UnNormalized.npy")
        norm_savename = os.path.join(self.encoding_output,
                                     f"{self.protein_name}_{self.encoding}_Normalized.npy")

        # Save the embeddings
        np.save(unnormalized_savename, unnormalized_embeddings)
        np.save(norm_savename, normalized_embeddings)
        
        
    # =========================================================================
    # ============== Protect instance variables as attributes =================
    # =========================================================================
    @property
    def encoding(self):
        return self._encoding

    @property
    def fasta_path(self):
        return self._fasta_path
        
    @property
    def target_protein_indices(self):
        return self._target_protein_indices

    @property
    def protein_name(self):
        return self._protein_name

    @property
    def output(self):
        return self._output

    @property
    def n_positions_combined(self):
        return self._n_positions_combined

    @property
    def combi_space(self):
        return self._combi_space

    @property
    def wt_seq(self):
        return self._wt_seq

    @property
    def wt_aas(self):
        return self._wt_aas

    @property
    def target_python_inds(self):
        return self._target_python_inds

    @property
    def encoding_output(self):
        return self._encoding_output

    @property
    def fasta_output(self):
        return self._fasta_output
    
    @property
    def all_combos(self):
        return self._all_combos
    
    @property
    def index_to_combo(self):
        return self._index_to_combo


if __name__=="__main__":
        
    # Import necessary modules and functions
    import argparse
    import os
    
    # Turn off extensive tensorflow readout
    os.environ['TF_CPP_MIN_LOG_LEVEL'] = '2' 
    
    # Import MLDE functions and classes
    from Support.Encode.SupportFuncs import calculate_batch_size, check_args

    # Instantiate argparser
    parser = argparse.ArgumentParser()

    # Add required arguments
    parser.add_argument("encoding", help = "Choice of 'onehot', 'georgiev', 'resnet', 'bepler', 'unirep', 'transformer', or 'lstm'")
    parser.add_argument("protein_name", help = "Protein name alias")
    parser.add_argument("--fasta", help = "FASTA file containing parent sequence", 
                        required = False, default = None, type = str)
#    parser.add_argument("--positions", help = "AA indices to target",
#                        required = False, nargs = "+", dest = "positions", default = None, type = str)
#    parser.add_argument("--n_combined", help = "Number of positions to combine",
#                        required = False, default = None, type = int)
    parser.add_argument("--output", help = "Save location for output files.",
                        required = False, default = os.getcwd())
#    parser.add_argument("--batches", help = "Number of batches for embedding calculations",
#                        required = False, type = int, default = 0)

    # Parse the arguments
    args = parser.parse_args()
    
    # Make sure the arguments are appropriate
    #check_args(args)    
 
    # Construct the embedding generator
    # without target_protein_indices nor n_positions_combined arguments
    embedding_obj = EncodingGenerator(args.encoding, args.protein_name, 
                                      fasta_path = args.fasta,
                                      target_protein_indices = None,
                                      n_positions_combined = None,
                                      output = args.output)

