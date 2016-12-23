from scipy import stats

# f = open("EMLIKE.txt", "r")
# l = f.readlines()
index = 0
l1 = range(600)
with open("likeBlocks.txt") as infile:
    for line1 in infile:
        l2 = map(lambda x: float(x), str.split(line1))
        slope, intercept, r_value, p_value, std_err = stats.linregress(l1, l2)
        # print slope, l[index]
        print slope
        index = index + 1
