"""
Модуль анализа последовательностей LCR
"""

from collections import Counter

# считаем представленность аминокислот
lcr = "GVGGVPGEVDGVTLGVVGVKMGVGGGGELKVVGLGGIISGGGGGSEAGADGGEDEDAD"
lcr2 = "GVGGVPGEVDGMGVGGGGESGGGGGSEAGADGGEDEDAD"
lcr3 = "GVSRSRSR"

def aa_counter(row):
    len_str = len(row)
    aa_lst = list(row)
    c = Counter(aa_lst).most_common(2)
    if (c[0][1]/len_str)/(c[1][1]/len_str) >=2.5:
        print(c, c[0][1]/len_str, c[1][1]/len_str, (c[0][1]/len_str)/(c[1][1]/len_str))
        return c[0][0]
    return None


def diaa_counter(row):
    """
    frequency of the second most abundant amino acid being more than half the frequency of the most abundant amino acid
    """

    len_str = len(row)
    aa_lst = list(row)
    c = Counter(aa_lst).most_common(2)
    if c[1][1]/len_str > (c[0][1]/len_str)/2:
        print(c, c[0][1]/len_str, c[1][1]/len_str, (c[0][1]/len_str)/(c[1][1]/len_str))
        return c[0][0]+c[1][0]


print(aa_counter(lcr))
print(aa_counter(lcr2))
print(diaa_counter(lcr3))