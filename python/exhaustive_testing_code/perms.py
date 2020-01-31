from math import factorial
from functools import reduce
from operator import mul
from itertools import accumulate as PrefixSum #prefix sum
from time import process_time

import cProfile


# This is proof of concept code.

# It evolved from a simple fixed script to one that exhaustively tests
# combinations of buckets and bucket sizes by introducing global
# variables and probably other abuses. It's not good practice in
# general but this is not going to be expanded further (that will be
# done in C) so it's OK here even if pylint rightly throws wobblies.

# In a few places I've ignored pylint's standards for improved
# readability.



bkts = -1	 # num of buckets per row
BitsPerItem = -1 # each bucket contains ipb items (see below)

# collection of globals. Some of these are aliases, but those should
# be removed xxx
ipb = -1
BaseMask = -1
lce = -1
cdims = -1
BitStrLen = -1
CubeLUT = [-99]
perms = None

INTLERR = "Internal Error:"

CPRINTON = 0



def cprint(*args):
    "Conditional print"
    if CPRINTON:
        print(*args, sep=" ", end="\n")


def assert_between(val, lo, hi, FuncName):
    assert lo <= val <= hi, \
        "{0} val outside lo or hi in routine '{1}'. " \
        "val / lo / hi are, respectively: {2} / {3} / {4}" \
	        .format(INTLERR, FuncName, val, lo, hi)


def CalcDependentData():
    """
    Updates the global state of the rest of the global variables based
    on the values of 'bkts' and 'BitsPerItem'
    """

    global ipb, BaseMask, Ice, cdims, BitStrLen
    global lce

    assert_between(bkts, 2, 8, "CalcDependentData")
    assert_between(BitsPerItem, 2, 8, "CalcDependentData")

    ipb = 2 ** BitsPerItem	# items per bucket; 0..ipb-1
    BaseMask = ipb - 1	# 2^3 = 8 -> mask of ..00111

    # handy equivalences.
    lce = ipb			# length of cube's edge
    cdims = bkts		# dimension of cube
    BitStrLen = lce ** bkts 	# cube volume / entire length of bit string



#what if #BitsPerItem is less than buckets? xxx


def RangeInc(start: int, end: int):
    """
    Convenience function. Behaves exactly as normal range(a, b) but b
    is inclusive.
    """
    assert start <= end, INTLERR + "starts beyond end in RangeInc"
    return range(start, end + 1)



def extract_bit_fields(i: int, BitFieldWidth: int, HowMany: int, allbits=True):

    """For 'i', extracts 'HowMany' bitfields of 'BitFieldWidth'.

    The bitfields are assumed to be from the right and contiguous,
    so for i = 0bzzzzz,11,10,01,00 (commas for clarity) then
    extract_bit_fields(i, 2, 4) will produce [0b11, 0b10, 0b01, 0b00].

    If 'allbits' = True then the rightmost reminder (here zzzzz) must
    be allbits zeros and an error will be thrown if they aren't.
    If 'allbits' = False no such check is made.
    """
    assert (i >= 0) & (BitFieldWidth > 0) & (HowMany > 0), \
    	("{0} Bad args for extract_bit_fields: i / BitFieldWidth / HowMany " + \
            "respectively {1} / {2} / {3}") .\
            format(INTLERR, i, BitFieldWidth, HowMany)

    itmp = i
    mask = (2 ** BitFieldWidth) - 1
    hm = HowMany
    res = []
    while hm > 0:
        res = [itmp & mask] + res
        itmp >>= BitFieldWidth
        hm -= 1

    if allbits & itmp:
        ErrMsg = "{0} extract_bit_fields: 'allbits' was specified but " + \
                 "remnant of 'i' was nonzero (was {1}".format(INTLERR, itmp)
        raise Exception(ErrMsg)
    else:
        return res 





def ExtractIndexesOf(val: int):
    """Gets the indexes of where 'val' is within the cube"""
    return extract_bit_fields(val, BitsPerItem, bkts)



def RisingFactorial(i: int, n: int):
    """the first 'n' terms of factorial('i') but ascending
    eg. fact(5, 3) = 5 * 6 * 7
    <https://en.wikipedia.org/wiki/Falling_and_rising_factorials>"""
    return reduce(mul, range(i, i + n))



def FullCubePopcount(xcdims, xlce):
    return RisingFactorial(xlce, xcdims) // factorial(xcdims)




# Manually build permutation table, so we can compare it to my
# calculations xxx



def PopulatePermsDict():
    """
    Fills in the permutations dictionary. This is the brute force
    creation of the mapping from unencoded to encoded values. It's
    used only to verify my proposed method works.

    """
    global perms
    nitp = ipb ** bkts	# num items to permute
    nnatp = 0		# next number allocated to permutation
    perms = dict()
    for x in range(nitp):
        tmp = ExtractIndexesOf(x)
        ToAdd = tuple(sorted(tmp))
        if ToAdd not in perms:
            perms[ToAdd] = nnatp
            nnatp += 1

    assert len(perms) == FullCubePopcount(bkts, ipb), \
        "internal error: slow and quick perm counts don't match. Are " + \
        str(FullCubePopcount(bkts, ipb)) + " and " + str(len(perms))



def CalcCubeLUT():
    "calculates the cube lookup table"
    global CubeLUT, lce
    CubeLUT = []
    CubeLUT.append(list(RangeInc(0, lce)))  # This line...
    # ... should be unnecessary but is a workaround for a bug I think
    # is caused by RisingFactorial not defining 0! = 1. It's OK for now.
    for cdim in RangeInc(1, cdims):
        y = [0]
        for xlce in range(lce, 0, -1):
            y.append( FullCubePopcount(cdim, xlce) )
        CubeLUT.append(list(PrefixSum(y)))



def ixOf(dim: int, frust: int):
    """
    Indexes into the CubeLUT.  'dim' is the dimension, 'frust' is the
    index of the running total of the frustrum. See related
    LibreOffice documentation.
    frust may be from -1 to 'lce'.
    Value at -1 is always zero, so the first layer may always subtract
    0. That adds an lookup for this special case but removes any
    conditional which will be worth it later (branches are
    expensive).
    This can be fudged in C for free by adjusting the base of the array.
    """
    assert_between(dim, 0, cdims, "ixOf (testing cdims)")  # xxx why not cdims - 1 - explain!
    assert_between(frust, -1, lce - 1, "ixOf (testing frust)")
    cprint("orig dim:", dim, ", cube slice:", CubeLUT[dim-1])
    return CubeLUT[dim-1][frust+1]



def GetRankOverall(i: int):
    """ Get the number of 1-bits at or below 'i' """
    assert i < BitStrLen, "{0} i >= BitStrLen in GetRankOverall: " + \
    	"i is {1}, BitStrLen is {2}".format(INTLERR, i, BitStrLen)

    uixs = ExtractIndexesOf(i)  # uixs = unmodified indexes
    uixs.sort()

    mixs = [i - 1 for i in uixs]  # mixs = modified indexes
    mixs[len(mixs) - 1] += 1  #undo decrement of last index
    cprint("mixs: ", mixs)

    # ofs = offsets
    ofs = [-1] + mixs[:-1]  # [4, 5, 6] -> [-1, 4, 5]
    ofs[len(ofs) - 1] += 1	# undo decrement of last offset;
    				# [-1, 4, 5] ->  [-1, 4, 6]
    cprint("ofs: ", ofs)
    assert len(mixs) == len(ofs), \
    	"{0} len(mixs) &  len(ofs) are unequal; {1} / {2}" \
        			.format(INTLERR, mixs, ofs)

    ixsofs = list(zip(mixs, ofs))
    cprint()
    cprint(" i at root: ", i, ",Unmodified indexes at root: ", uixs)
    cprint(" ixsofs at root: ", ixsofs, "\n")
    return GetRankWithinCube(len(ixsofs), ixsofs)




def GetRankWithinCube(dim: int, ixofs):
    """
    Gets the rank of the 1 bit in the n-cube, that rank being the
    perm-packed value we want.  See associated LibreOffice docs which
    explain this clearly
    """
    assert dim == len(ixofs), "dim / len(ixofs) don't match: are {0} / {1}" \
        				.format(dim, len(ixofs))

    # This will work also if it is commented back in, and the 'if dim
    # == 1 ... return 0' below is commented out. It's a degenerate
    # case and this is just an optimisation, but a good one.

    # if dim == 0:
    #    return 0

    (ixof, *rest) = ixofs
    (frustix, offsetix) = ixof
    if dim == 1:
        return frustix - offsetix
    else:
        FrustTopSize = ixOf(dim, frustix)
        OffsetSize = ixOf(dim, offsetix)
        cprint("dim = {0}, (frustix, offsetix) = {1}, rest = {2}".\
               format(dim, ixof, rest))
        cprint("FrustTopSize ={0}, OffsetSize = {1}".\
               format(FrustTopSize, OffsetSize))
        FrustSize = FrustTopSize - OffsetSize
        cprint("FrustSize = {0}\n".format(FrustSize))

        SliceSize = GetRankWithinCube(dim - 1, rest)
        tot = FrustSize + SliceSize
        cprint("FrustSize, SliceSize, tot: {0}, {1}, {2}\n".\
              format(FrustSize, SliceSize, tot))

        return tot




def validate():
    """
    Checks by exhaustion that my method and the simple, brute-force
    method produce the same results.
    """
    for ctr in range(BitStrLen):
        ro = GetRankOverall(ctr)
        ixs = ExtractIndexesOf(ctr)
        perm = perms[tuple(sorted(ixs))]
        cprint("result -- perm={0}, ro={1}, ctr={2}, sorted(ixs) = {3}\n\n".\
               format(perm, ro, ctr, sorted(ixs)))
        assert perm == ro, (INTLERR + \
                    "\nmismatch: perm={0}, ro={1}, ctr={2}, sorted(ixs) = {3}").\
                            format(perm, ro, ctr, sorted(ixs))




        
def TestAll():
    test_extract_bit_fields()

# TestAll() xxx sort this later

for bkts in RangeInc(2, 8):
    for BitsPerItem in RangeInc(2, 8):
        CalcDependentData()
#        print("x" * 4096)  # force a flush of stdout, else buffers too much
        print("\n\nbkts = {0}, BitsPerItem = {1}, BitStrLen = {2}"\
              .format(bkts, BitsPerItem, BitStrLen))
        # below, if 2**24 then ~ 1 hour runtime, if 2**28 then overnight
        if BitStrLen > 2**24:  # excessively large combo?
            print("skipping...")
            continue
        toverall = process_time()
        CalcCubeLUT()
        cprint(CubeLUT)
        cprint(CubeLUT, "\n")
        PopulatePermsDict()
        cprint(perms)
        tpopdict = process_time()
        print("    populated dict time: ", int(tpopdict - toverall))
        validate()
        tvalidate = process_time()
        print("    validate time: ", int(tvalidate - tpopdict))
        print("overall time: ", int(tvalidate - toverall))







def test_extract_bit_fields():
    # test for 2 bit fields, extract 4
    tst = 0b11100100
    expected2x4 = [0b11, 0b10, 0b01, 0b00]
    y = extract_bit_fields(tst, 2, 4)
    if y != expected2x4:
        raise Exception(INTLERR + " failed to extract bits correctly")

    tst = 0b111100100
    DidNotExcept = True
    try:
        y = extract_bit_fields(tst, 2, 4)
    except:
        DidNotExcept = False
    if DidNotExcept:
        raise Exception(INTLERR + "failed to except; remnant expected")


    # and test for 5 bit fields, extract 3
    tst = 0b110010101100110
    expected5x3 = [0b11001, 0b01011, 0b00110]
    y = extract_bit_fields(tst, 5, 3)
    if y != expected5x3:
        raise Exception(INTLERR + " failed to extract bits correctly")

    tst = 0b1110010101100110
    DidNotExcept = True
    try:
        y = extract_bit_fields(tst, 2, 4)
    except:
        DidNotExcept = False
    if DidNotExcept:
        raise Exception(INTLERR + "failed to except; remnant expected")

