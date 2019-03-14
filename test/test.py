# import seqtools
# from os import path
# x = (None, False, True)

# print(x)
# print(len(x))

# if x[0]:
#     print(f"x[0]: {x[0]}")

# if x[1]:
#     print(f"x[1]: {x[1]}")

# if x[2]:
#     print(f"x[2]: {x[2]}")

# bla = [0,0,1,1]

# print(len(bla), sum(bla))

# test = "axbabdfadfbdfagfdavdsavsdav"
# test1 = [0,1,2,3,4,5,6,7,8,9,10]

# for x,y in zip(test, test1):
#     print(x,y)


from seqtools.generate import generate_dna

generate_dna(30)
