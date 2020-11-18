from ase.collections import g2


def get_permutations(array):
    ''' Get all possible permutaiton for a give array '''
    valid_permut = []
    index1 = 0
    permut_helper(index1, array, valid_permut)
    valid_permut = [''.join(permut) for permut in valid_permut]
    return valid_permut


def permut_helper(index1, array, valid_permut):
    if index1 == len(array) - 1:
        # array[:] return all elements with sliced operation
        # and will get rid of repeated permutations
        valid_permut.append(array[:])
    for index2 in range(index1, len(array)):
        swap(array, index1, index2)
        permut_helper(index1 + 1, array, valid_permut)
        swap(array, index1, index2)


def swap(array, index1, index2):
    array[index1], array[index2] = array[index2], array[index1]
