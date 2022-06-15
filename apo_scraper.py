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


def check_related_pdb_ids_apo(pdb_id):
    """Obtains all related proteins for a given PDB and checks if they are apo

    Args:
        pdb_id (str): 4 letter ID for structure

    Returns:
        (dict): Dictionary, for the ID checked as key, with "apo" structures
        as values
    """
    apo_dict = {}
    all_info = get_all_pdb_info(pdb_id)
    if not all_info:
        # RAISE ERROR
        pass
    else:
        for i in range(len(all_info['pdbx_database_related'])):
            new_pdb = all_info['pdbx_database_related'][i]['db_id']
            is_apo = check_if_apo(new_pdb)
            if is_apo:
                append_value(apo_dict, pdb_id, new_pdb)
    return apo_dict


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb_id', type=str,
                        help='ID of PDB structure from RCSB PDB that you want'
                             'to find apo structures of')
    args = parser.parse_args()
    print(check_related_pdb_ids_apo(args.pdb_id))
