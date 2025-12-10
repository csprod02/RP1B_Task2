from unittest import TestCase, main
# Will start with just this one
from mutate_genome import section_overlapping, check_for_deletion, apply_mutations


## TEST CLASS 1
class Test_deletions_overlapping(TestCase):
    # mock deletion array to test cases against
    del_array = [{'start_index': 5, 'length': 4}, # deletion covers 5 -> 8 (inclusive)
                {'start_index': 28, 'length': 10}] # deletion covers 28 -> 37 (inclusive)

    def test_overlapping_CASE_1(self): # end of new deletion overlaps existing one
        new_start = 2
        new_length = 7 # deletion ends at index 6
        
        expected = True
        result = section_overlapping(new_start, new_length, self.del_array)
        
        self.assertEqual(result, expected)

    def test_overlapping_CASE_2(self): # start of new deletion overlaps existing one
        new_start = 7
        new_length = 6 # deletion ends at index 12
        
        expected = True
        result = section_overlapping(new_start, new_length, self.del_array)
        
        self.assertEqual(result, expected)

    def test_overlapping_CASE_3(self): # new deletion starts and ends within existing one
        new_start = 30
        new_length = 3 # deletion ends at index 32
        
        expected = True
        result = section_overlapping(new_start, new_length, self.del_array)
        
        self.assertEqual(result, expected)

    def test_overlapping_CASE_4(self): # new deletion starts and ends outside of existing deletion, but encompases it entirely
        new_start = 3
        new_length = 10 # deletion ends at index 12
        
        expected = True
        result = section_overlapping(new_start, new_length, self.del_array)
        
        self.assertEqual(result, expected)

    def test_overlapping_CASE_5(self): # no overlap found
        new_start = 10
        new_length = 6
        
        expected = False
        result = section_overlapping(new_start, new_length, self.del_array)
        
        self.assertEqual(result, expected)



# TEST CLASS 2
class Test_check_for_deletion(TestCase):
    del_array = [{'start_index': 5, 'length': 4}, # deletion covers 5 -> 8 (inclusive)
                {'start_index': 28, 'length': 10}] # deletion covers 28 -> 37 (inclusive)
    
    def test_CASE_1(self): # suggested SNP location is inside of a deletion
        base_pos = 6

        expected = True
        result = check_for_deletion(base_pos, self.del_array)

        self.assertEqual(result, expected)

    def test_CASE_2(self): # suggested SNP location is first index of an existing deletion
        base_pos = 5

        expected = True
        result = check_for_deletion(base_pos, self.del_array)

        self.assertEqual(result, expected)

    def test_CASE_3(self): # suggested SNP location is last index of an existing deletion
        base_pos = 37

        expected = True
        result = check_for_deletion(base_pos, self.del_array)

        self.assertEqual(result, expected)

    def test_CASE_4(self): # suggested SNP location is directly after an existing deletion
        base_pos = 38

        expected = False
        result = check_for_deletion(base_pos, self.del_array)

        self.assertEqual(result, expected)

    def test_CASE_5(self): # suggested SNP location is directly before an existing deletion
        base_pos = 27

        expected = False
        result = check_for_deletion(base_pos, self.del_array)

        self.assertEqual(result, expected)

    def test_CASE_6(self): # SNP location is not in or next to a deletion zone
        base_pos = 15
        
        expected = False
        result = check_for_deletion(base_pos, self.del_array)

        self.assertEqual(result, expected)




## TEST CLASS 3
class Test_apply_mutations(TestCase):
    test_genome = '0123456789'

    def test_CASE_1(self): # applying a single deletion
        mutations_list = [{'type': 'deletion',
                          'start_index': 2,
                          'length': 4}]

        expected = '016789'
        result = apply_mutations(self.test_genome, mutations_list)

        self.assertEqual(result, expected)

    def test_CASE_2(self): # applying a single SNP
        mutations_list = [{'type': 'snp',
                          'base_index': 4,
                          'alt': 'X'}]

        expected = '0123X56789'
        result = apply_mutations(self.test_genome, mutations_list)

        self.assertEqual(result, expected)
    
    def test_CASE_3(self): # applying a single insertion
        mutations_list = [{'type': 'insertion',
                          'insertion_string': 'XYZABC',
                          'insert_loc': 7}]

        expected = '01234567XYZABC89'
        result = apply_mutations(self.test_genome, mutations_list)

        self.assertEqual(result, expected)

    def test_CASE_4_1(self): # deletion -> snp -> insertion
        mutations_list = [{'type': 'deletion', # covers indexes 2345
                          'start_index': 2,
                          'length': 4},
                         
                         {'type': 'snp',
                          'base_index': 7,
                          'alt': 'X'},
                         
                         {'type': 'insertion', # inserted between indexes 8 and 9
                          'insertion_string': 'XYZABC',
                          'insert_loc': 8}]

        sorted_mutations = sorted(mutations_list, key= lambda x : x.get('start_index', x.get('base_index', x.get('insert_loc'))), reverse=True)
        
        expected = '016X8XYZABC9'
        result = apply_mutations(self.test_genome, sorted_mutations)

        self.assertEqual(result, expected)

    def test_CASE_4_2(self): # deletion -> insertion -> snp
        mutations_list = [{'type': 'deletion', # covers indexes 2345
                          'start_index': 2,
                          'length': 4},
                         
                         {'type': 'insertion', # inserted between indexes 6 and 7
                          'insertion_string': 'XYZABC',
                          'insert_loc': 6},
                         
                         {'type': 'snp',
                          'base_index': 7, # first index after insertion should be modified
                          'alt': 'X'}]

        sorted_mutations = sorted(mutations_list, key= lambda x : x.get('start_index', x.get('base_index', x.get('insert_loc'))), reverse=True)
        
        expected = '016XYZABCX89'
        result = apply_mutations(self.test_genome, sorted_mutations)

        self.assertEqual(result, expected)

    def test_CASE_4_3(self): # insertion -> snp -> deletion
        mutations_list = [{'type': 'insertion', # inserted between indexes 3 and 4
                          'insertion_string': 'XYZABC',
                          'insert_loc': 3},
                         
                         {'type': 'snp',
                          'base_index': 6,
                          'alt': 'X'},
        
                         {'type': 'deletion', # covers indexes 2345
                          'start_index': 8,
                          'length': 2}]

        sorted_mutations = sorted(mutations_list, key= lambda x : x.get('start_index', x.get('base_index', x.get('insert_loc'))), reverse=True)
        
        expected = '0123XYZABC45X7'
        result = apply_mutations(self.test_genome, sorted_mutations)

        self.assertEqual(result, expected)

    def test_CASE_4_4(self): # insertion -> deletion -> snp
        mutations_list = [{'type': 'insertion', # inserted between indexes 3 and 4
                          'insertion_string': 'XYZABC',
                          'insert_loc': 3},
        
                         {'type': 'deletion', # covers indexes 4567
                          'start_index': 4,
                          'length': 4},
                         
                         {'type': 'snp',
                          'base_index': 9,
                          'alt': 'X'},]

        sorted_mutations = sorted(mutations_list, key= lambda x : x.get('start_index', x.get('base_index', x.get('insert_loc'))), reverse=True)
        
        expected = '0123XYZABC8X'
        result = apply_mutations(self.test_genome, sorted_mutations)

        self.assertEqual(result, expected)

    def test_CASE_4_5(self): # snp -> deletion -> insertion
        mutations_list = [{'type': 'snp',
                          'base_index': 0,
                          'alt': 'X'},

                         {'type': 'deletion', # covers indexes 4567
                          'start_index': 4,
                          'length': 4},
                          
                         {'type': 'insertion', # inserted after index 9
                          'insertion_string': 'XYZABC',
                          'insert_loc': 9}]

        sorted_mutations = sorted(mutations_list, key= lambda x : x.get('start_index', x.get('base_index', x.get('insert_loc'))), reverse=True)
        
        expected = 'X12389XYZABC'
        result = apply_mutations(self.test_genome, sorted_mutations)

        self.assertEqual(result, expected)

    def test_CASE_4_6(self): # snp -> insertion -> deletion
        mutations_list = [{'type': 'snp',
                          'base_index': 0,
                          'alt': 'X'},
                          
                         {'type': 'insertion', # inserted between indexes 5 and 6
                          'insertion_string': 'XYZABC',
                          'insert_loc': 5},
                         
                         {'type': 'deletion', # covers indexes 78
                          'start_index': 7,
                          'length': 2}]

        sorted_mutations = sorted(mutations_list, key= lambda x : x.get('start_index', x.get('base_index', x.get('insert_loc'))), reverse=True)
        
        expected = 'X12345XYZABC69'
        result = apply_mutations(self.test_genome, sorted_mutations)

        self.assertEqual(result, expected)



## MAIN
if __name__ == '__main__':
    main()
    
