from libtbx import easy_pickle
prefix="/Users/wyf/Desktop/align_code/myDB/antinon3"
codes=easy_pickle.load(prefix+".codes")
nlm_coefs=easy_pickle.load(prefix+".nlm")
nn_coefs=easy_pickle.load(prefix+".nn")
nl_coefs=easy_pickle.load(prefix+".nl")
rmaxs=easy_pickle.load(prefix+".rmax")
# rmaxs=easy_pickle.load("~/database")

#retrieval information for a specific model
# code="CHEMBL521242" 
# indx=codes.index(code)
# for name in codes:
# 	print codes.index(name)
indx=25121
code=codes[indx]
# nlm_coef=nlm_coefs[indx]
# nn_coef=nn_coefs[indx]
# nl_coef=nl_coefs[indx]
rmax=rmaxs[indx]

# for i in range(len(nl_coef)):
# 	print "nl_coef[", i, "]: ",nl_coef[i]
# for j in range(len(nlm_coef)):
# 	print "nlm_coef[", j,"]: ",nlm_coef[i]
# print "--------------------------------------------------"
# print "nl_coef[", i, "]  : ",nl_coef[i]
# print "nlm_coef[", j,"]: ",nlm_coef[j]
# print "codes :      "
# print "len(nl_coefs) :      ", len(nl_coefs)
# print "len(nlm_coefs):      ", len(nlm_coefs)
print"indx:       ", indx
print"code:       ", code
print"rmax:       ", rmax