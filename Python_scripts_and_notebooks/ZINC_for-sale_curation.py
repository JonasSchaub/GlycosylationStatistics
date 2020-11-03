"""
TODO
http://www.rdkit.org/docs/GettingStartedInPython.html#picking-diverse-molecules-using-fingerprints
"""

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker


class ZincPicker:
    """
    TODO
    500 iterations where 1000 molecules are picked out of 44,000 each iteration makes 500,000 molecules picked out of
    22,000,000 in total.
    """
    # static variables
    pool_size = 44000
    pick_size = 1000
    num_iterations = 500
    input_file_path = 'ZINC_for-sale_subset_2020_Okt_19.txt' # 22,252,431 molecules in this set
    output_file_path = 'output/ZINC_for-sale_picked_subset.txt'

    def __init__(self):
        """
        TODO
        """
        # instance variables
        self.fingerprint_list = []

    def calculate_dice_similarity_distance(self, i, j):
        """
        TODO
        Function to calculate the distance between two molecular fingerprints from a list using dice similarity.

        :param i:
        :param j:
        :param fps:
        :return:
        :rtype: object
        """
        return 1 - DataStructs.DiceSimilarity(self.fingerprint_list[i], self.fingerprint_list[j])


    def execute(self):
        """
        TODO
        """
        print()
        print("Loading input file with path: " + ZincPicker.input_file_path)
        zinc_for_sale_mol_supplier = Chem.SmilesMolSupplier(ZincPicker.input_file_path)
        num_none_mols = 0
        print("Output file path: " + ZincPicker.output_file_path)
        writer = Chem.SmilesWriter(ZincPicker.output_file_path)
        lower_index = 0
        upper_index = ZincPicker.pool_size
        print("Entering picking iterations...")
        print()
        for y in range(0, ZincPicker.num_iterations):
            print("Number of iteration: ", y)
            print("Lower index: ", lower_index)
            print("Upper index: ", upper_index)
            print("Loading molecules now...")
            molecules = []
            for x in range(lower_index, upper_index):
                mol = zinc_for_sale_mol_supplier[x]
                if mol is None:
                    num_none_mols += 1
                    continue
                molecules.append(mol)
            while molecules.count(None):
                molecules.remove(None)
            # radius 3
            print("Number of molecules loaded: ", len(molecules))
            print("Calculating fingerprints...")
            self.fingerprint_list = [GetMorganFingerprint(x, 3) for x in molecules]
            nfps = len(self.fingerprint_list)
            print("Number of fingerprints: ", nfps)
            print("Now min-max picking ", ZincPicker.pick_size, " out of the finger print list...")
            picker = MaxMinPicker()
            pickIndices = picker.LazyPick(self.calculate_dice_similarity_distance, nfps, ZincPicker.pick_size, seed=23)
            print("Finished picking, writing to file...")
            for z in pickIndices:
                writer.write(molecules[z])
            # clear memory
            molecules = []
            self.fingerprint_list = []
            nfps = 0
            picker = None
            pickIndices = []
            # raise indices
            lower_index = lower_index + ZincPicker.pool_size
            upper_index = upper_index + ZincPicker.pool_size
            print("Finished this iteration, entering the next...")
            print()
        print("Execution successful.")
        print("Picked ", ZincPicker.pick_size * ZincPicker.num_iterations - num_none_mols, " out of ",
              ZincPicker.num_iterations * ZincPicker.pool_size, " molecules in ", ZincPicker.num_iterations,
              " iterations, while picking ", ZincPicker.pick_size, " in each iteration.")

if __name__ == '__main__':
    picker_script = ZincPicker()
    picker_script.execute()
