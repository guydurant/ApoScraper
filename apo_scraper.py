from pypdb.clients.search.operators.sequence_operators import SequenceOperator
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.fasta.fasta_client import get_fasta_from_rcsb_entry
from pypdb import get_info, get_pdb_file
from non_ligand_chemicals import NON_LIGAND_HETATMS
import warnings
import argparse
warnings.filterwarnings("ignore", category=DeprecationWarning)


def get_all_pdb_info(pdb_id):
    """Gets all information about a given PDB structures

    Args:
        pdb_id (str): 4 letter ID for structure

    Returns:
        (str): All information about the PDB structure using pypdb package
    """
    return get_info(pdb_id)


def check_if_apo(pdb_id):
    """Checks if given PDB code relates to an apo structure

    Args:
        pdb_id (str): 4 letter ID for structure

    Returns:
        (bool): True if apo structure, False if not
    """
    pdb_file = get_pdb_file(pdb_id, filetype='pdb', compression=False)
    # Converts all parts of the file into a list and removes spaces
    cleaned_pdb_file = list(filter(('').__ne__, pdb_file.split(' ')))
    # Gets indices for all places heteroatoms are mentioned
    raw_coord = [i for i, j in enumerate(
        cleaned_pdb_file) if j == "\nHETNAM" or j == "\nHETSYN"]
    ligs = []
    # Iterates through the indices (along one) to get 3 letter code for
    # chemicals
    for i in [a + 1 for a in raw_coord]:
        # Checks if 3 letter code not a NON_LIGAND_HETATMS
        if cleaned_pdb_file[i] not in NON_LIGAND_HETATMS:
            # Sometimes weird strings are in wrong position so have to check
            # one along
            if len(cleaned_pdb_file[i]) < 3:
                if cleaned_pdb_file[i+1] not in NON_LIGAND_HETATMS:
                    ligs.append(cleaned_pdb_file[i+1])
            else:
                ligs.append(cleaned_pdb_file[i])
    # Stores ligands in list, currently pointless but useful in future
    # If ligs has nothing in it, must be an apo structure, else it is not
    if len(ligs) > 0:
        return False
    else:
        return True


def append_value(dict_obj, key, value):
    """Helper function to append to dictionaries

    Check if key exist in dict or not, if it does, makes value a list and
    appends to it. If not, adds to dictionary normally.

    Args:
        dict_obj (dict): dictoinary appending into
        key (any): key
        value (any): value
    """
    if key in dict_obj:
        # Key exist in dict.
        # Check if type of value of key is list or not
        if not isinstance(dict_obj[key], list):
            # If type is not list then make it list
            dict_obj[key] = [dict_obj[key]]
        # Append the value in list
        dict_obj[key].append(value)
    else:
        # As key is not in dict,
        # so, add key-value pair
        dict_obj[key] = value


def check_related_pdb_ids_apo(pdb_id, seq):
    """Obtains all related proteins for a given PDB and checks if they are apo

    Args:
        pdb_id (str): 4 letter ID for structure
        seq (float): Sequence similarity to search similar structures by.
        Between 0 and 1.

    Returns:
        (dict): Dictionary, for the ID checked as key, with "apo" structures
        as values
    """
    apo_dict = {}
    similar_pdbs = sequence_similarity_search_for_pdb_id(pdb_id, seq)
    for pdb in similar_pdbs:
        try:
            is_apo = check_if_apo(pdb)
        except AttributeError:
            print(f'{pdb} check failed')
        if is_apo:
            append_value(apo_dict, pdb_id, pdb)
    return apo_dict


def sequence_similarity_search_for_pdb_id(pdb_id, seq):
    """Performs a sequence similairty search for a given PDB structure and a
    given threshold

    Args:
        pdb_id (str): 4 letter ID for structure
        seq (float): Sequence similarity to search similar structures by.
        Between 0 and 1.

    Returns:
        list: All similiar PDBS to given PDB ID
    """
    # Fetches the first sequence in the PDB ID fasta file
    fasta_sequence = get_fasta_from_rcsb_entry(pdb_id)[0].sequence

    # Performs sequence search ('BLAST'-like) using the FASTA sequence
    results = perform_search(
        return_type=ReturnType.ENTRY,
        search_operator=SequenceOperator(
            sequence=fasta_sequence,
            identity_cutoff=seq,
            evalue_cutoff=1000
        ),
        return_with_scores=True
    )
    return process_sequence_similarity_search_results(results)


def process_sequence_similarity_search_results(results):
    """Converts the results of the sequence search into a list of PDB ids

    Args:
        results (list): List of ScoredResult classes

    Returns:
        list: List of PDB ids
    """
    if len(results) > 0:
        pdb_ids = [i.entity_id for i in results]
    return pdb_ids


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_id', type=str,
                        help='ID of PDB structure from RCSB PDB that you want'
                             'to find apo structures of')
    parser.add_argument('seq', type=float,
                        help='Sequence similarity threshold to search apo'
                             'structures for. Value must be between 0 and 1.')
    args = parser.parse_args()
    print(check_related_pdb_ids_apo(args.pdb_id, args.seq))
