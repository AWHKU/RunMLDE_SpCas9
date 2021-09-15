"""
This module contains a class needed to generate all encoding types.
"""
# Import necessary MLDE modules
from Support.Encode.MolBioInfo import all_aas, allowed_aas
from Support.Encode.GeorgievParams import georgiev_parameters
from Support.Encode.TapeModelLocations import tape_model_locations

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
print(os.path.dirname(os.path.realpath('__file__')))
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
    flat_encodings = np.reshape(unnormalized_encodings,[len(unnormalized_encodings), flat_shape])
    # Mean-center and unit variance scale the embeddings.
    means = flat_encodings.mean(axis=0)
    stds = flat_encodings.std(axis=0)
    normalized_flat_encodings = (flat_encodings - means)/stds    
    # Reshape the normalized embeddings back to the original form
    normalized_encodings = np.reshape(normalized_flat_encodings,unnormalized_encodings.shape)    
    return normalized_encodings

# Define a class for generating embeddings
class EncodingGenerator():
    """
    Directly input a set of fasta sequences from the mutagensis screen and send to make learnt embedding. 
    """
    # Initialize the embedding
    def __init__(self, encoding, protein_name,
                 fasta_path = None, target_protein_indices = None, 
                 n_positions_combined = None, output = os.getcwd()):
        # Assign all inputs as instance variables
        self._encoding = encoding.lower()
        self._fasta_path = fasta_path
        self._target_protein_indices = None
        self._protein_name = protein_name
        self._output = output      
        # Additional checks if working with learned embeddings. 
        if encoding in {"resnet", "bepler", "unirep", "transformer", "lstm"}:
            
            # Assert that we have the correct variables present
#            assert target_protein_indices is not None, "Did not define target indices"
            assert fasta_path is not None, "Did not specify location of fasta file"            
#            # Load the fasta file
#            self._process_input_fasta()
#
#            # Check the input indices
#            self._check_indices()        
#        elif encoding in {"georgiev", "onehot"}:
#            
#            # Assert that we have the correct variables present
#            assert n_positions_combined is not None, "Did not define n_positions_combined"
#            
#            # Assign a combinatorial space size and get a count of the variant size
#            self._n_positions_combined = n_positions_combined        
        else:
            raise AssertionError("Unknown encoding")        
        # Define the size of the combinatorial space
#        self._combi_space = 20**self.n_positions_combined        
        # Build output directories
        self._build_output_dirs()       
        # Build the list of combinations for the position and the dictionaries 
        # linking position index to combo
        #self._build_combo_dicts()
        self._generate_encodings()
    # Write a function that builds output directories
    def _build_output_dirs(self):
        """
        Self-explanatory: Build necessary directories for saving data.
        """
        # Get the start time
        init_time = strftime("%Y%m%d-%H%M%S")       
        # Build the output directories consistent for all encodings
        self._encoding_output = os.path.join(self._output, init_time, "Encodings")
        os.makedirs(self._encoding_output)        
        # Build the output directories only used only for generating learned embeddings
        self._fasta_output = os.path.join(self._output, init_time, "Fastas")
        os.makedirs(self._fasta_output)
            
   # Write a function that produces dictionaries linking combo and index in the
    # output encodings
    def _build_combo_dicts(self):
        """
        Builds dictionaries which link the identity of a combination (e.g. ACTV)
        to that combination's index in an encoding array, and vice versa. Both
        dictionaries are saved to disk.
        Now we rewrite this function to take in the Combo_setting from the fasta file, according to the order in the fasta file.
        """
        # Identify all possible combinations
        #self._all_combos = list(product(all_aas, repeat = self.n_positions_combined))        
        # Link combo to index and index to combo
        combo_to_index = {"".join(combo): i for i, combo in enumerate(self._all_combos)}
        self._index_to_combo = {i: "".join(combo) for i, combo in enumerate(self._all_combos)}        
        # Save the constructed dictionaries
        with open(os.path.join(self._encoding_output, f"{self._protein_name}_{self._encoding}_ComboToIndex.pkl"), "wb") as f:
            pickle.dump(combo_to_index, f)
        with open(os.path.join(self._encoding_output, f"{self._protein_name}_{self._encoding}_IndexToCombo.pkl"), "wb") as f:
            pickle.dump(self._index_to_combo, f)
                
    # Write a function that generates fasta files for the requested variants
    def _build_fastas(self):
        """
        Modify code to split input fasta sequences into batches of 20 sequences? 
        Then feed the sequences to TAPE
        """
        full_fasta=[]
        full_fasta_header=[]
        with open(self._fasta_path, "r+") as f:
            data = f.readlines()
            for line1, line2 in list(zip(data, data[1:]))[::2]:
                full_fasta.append("".join([line1, line2,]))
                full_fasta_header.append(line1)
        ### build WTaa profile
        WTaa=full_fasta_header[0].split("_")
        WTaa=WTaa[1].strip()
        WTaa=WTaa.split("-")
        WTaa=[aa[:-1] for aa in WTaa]
        ### make self._all_combos here
        combo_list=[]
        for j in full_fasta_header:
            j=j.strip()
            combos=j.split("_")[-1]
            combos=combos.split("-")
            combos=[m[-1] for m in combos]
            combo_list.append(combos)
        self._all_combos = combo_list
        self._build_combo_dicts()
        # Create a list to store all fasta filenames in
        n_batches = len(full_fasta)//20 
        if len(full_fasta)%20 >0:
            n_batches+=1
        fasta_filenames = [None for _ in range(n_batches)]        
        fasta_counter=0
        for i in range(n_batches):
            # Create a filename for the file we will use to store fasta data
            fasta_filename = os.path.join(self._fasta_output, "{}_Batch{}_Variants.fasta".format(self._protein_name, i))
            # Record the fasta_filename
            fasta_filenames[i] = fasta_filename            
            # Create a list in which we will store SeqRecords
            temp_seqs = full_fasta[fasta_counter:fasta_counter+20]
            fasta_counter+=20            
            # Write fasta sequences of the variant combinations to the file
            with open(fasta_filename, "w") as f:
                f.write("".join(temp_seqs))
        # Return the filename of the fasta sequences
        return fasta_filenames, n_batches, combo_list, WTaa
     
    # Write a function that generates learned encodings
    def _generate_learned(self):
        """
        Encodes a given combinatorial space using tape.
        Unlike Georgiev and one-hot encodings, these encodings are context-
        aware. To save on system RAM, this task is split into n_batches number
        of batches.
        """
        # Build fasta files
        fasta_filenames, n_batches, combo_list, WTaa = self._build_fastas()
        # make index for mutation locations
        self.target_protein_indices=WTaa
        index_checks = [_aa_ind_splitter.match(index) for index in self.target_protein_indices]
        aa_letter_index = [[match.group(1), int(match.group(2)) - 1] for match in index_checks]
        self._target_python_inds = [ind for _, ind in aa_letter_index]
        self._all_combos = combo_list
        # Create a list to store the names of the raw embedding files
        extracted_embeddings = [None for _ in range(n_batches)]
        # Get the name of the model and the load-from name
        weight_loc = tape_model_locations[self._encoding]       
        # Get a temporary filename for storing batches
        temp_filename = os.path.join(_filedir, "TempOutputs.pkl")
        # Loop over the number of batches
        for i, fasta_filename in tqdm(enumerate(fasta_filenames),desc = "Batch#", total = len(fasta_filenames),position = 0):
            # Run TAPE to get the transformer embeddings
            _ = subprocess.run(["tape-embed", fasta_filename, self._encoding, "--load-from", weight_loc, "--output", temp_filename])
            # Load the raw embedding file that was generated
            with open(temp_filename, "rb") as f:
                 raw_embeddings = pickle.load(f)            
            # Extract just the indices we care about
            extracted_embeddings[i] = np.array([protein_embedding[0, self._target_python_inds, :] for protein_embedding in raw_embeddings])       
        # Delete the temporary outputs file
        os.remove(temp_filename)        
        # Return the extracted embeddings, concatenating along the first dimension
        return np.concatenate(extracted_embeddings)
               
    # Write a function that generates encodings
    def _generate_encodings(self):
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
        if self._encoding in {"resnet", "bepler", "unirep", "transformer", "lstm"}:unnormalized_embeddings = self._generate_learned()
        else:
            raise ValueError(f"Unknown encoding type {self._encoding}")        
        # Normalize embeddings
        # Reshape the normalized embeddings back to the original form
        normalized_embeddings = _normalize_encodings(unnormalized_embeddings)
        # Create filenames for saving the embeddings
        unnormalized_savename = os.path.join(self._encoding_output, f"{self._protein_name}_{self._encoding}_UnNormalized.npy")
        norm_savename = os.path.join(self._encoding_output, f"{self._protein_name}_{self._encoding}_Normalized.npy")
        # Save the embeddings
        np.save(unnormalized_savename, unnormalized_embeddings)
        np.save(norm_savename, normalized_embeddings)                



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


