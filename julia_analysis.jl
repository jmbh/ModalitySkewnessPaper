using CSV
using MixedModels
using DataFrames
using CategoricalArrays

df = CSV.read("Files/ModalityML_Data.csv", DataFrame)
df.study = categorical(df.study)
df.par = categorical(df.par)

##################################################################
### Start with a base model and make it more complex step by step
##################################################################

# Base model
form0 = @formula(M ~ 1 + (1 | study) + (1 | par))

# Level-1 Predictors
form1 = @formula(M ~ 1 + valence + (1 | study) + (1 | par))
form1a = @formula(M ~ 1 + valence + (1 + valence | study) + (1 + valence | par))

# Level-2 Predictors
form2 = @formula(M ~ 1 + valence + TSlength_c + (1 + valence | study) + (1 + valence | par))
form2a = @formula(M ~ 1 + valence * TSlength_c + (1 + valence | study) + (1 + valence | par))
#form2aa = @formula(M ~ 1 + valence * TSlength_c + (1 | study) + (1 + valence | par))
#form2aaa = @formula(M ~ 1 + valence * TSlength_c + (1 + valence | study) + (1 | par))

form2b = @formula(M ~ 1 + valence * TSlength_c + zerocorr(1 + valence | study) + (1 + valence | par))
form2c = @formula(M ~ 1 + valence * TSlength_c + (1 + valence + TSlength_c | study) + (1 + valence + TSlength_c | par))


# Level-3 Predictors
form3 = @formula(M ~ 1 + valence * TSlength_c + retro + (1 + valence | study) + (1 + valence | par))
form3a = @formula(M ~ 1 + valence * TSlength_c + retro + zerocorr(1 + valence | study) + (1 + valence | par))
form3b = @formula(M ~ 1 + valence * TSlength_c + valence * retro + (1 + valence | study) + (1 + valence | par))

form4 = @formula(M ~ 1 + valence * TSlength_c + scale + (1 + valence | study) + (1 + valence | par))
form4a = @formula(M ~ 1 + valence * TSlength_c + valence * scale + (1 | study) + (1 + valence | par))

form4b = @formula(M ~ 1 + valence * TSlength_c + scale + zerocorr(1 + valence | study) + (1 + valence | par))

form5 = @formula(M ~ 1 + valence * TSlength_c + MpD_c + (1 + valence | study) + (1 + valence | par))
form5a = @formula(M ~ 1 + valence * TSlength_c + valence * MpD_c + (1 + valence | study) + (1 + valence | par))

form6 = @formula(M ~ 1 + valence * TSlength_c + valence * retro + valence * scale + valence * MpD_c + (1 | study) + (1 + valence | par))



m0 = fit(MixedModel, form0, df, Bernoulli())
m1 = fit(MixedModel, form1, df, Bernoulli())
m1a = fit(MixedModel, form1a, df, Bernoulli())

bic(m0)
bic(m1)
bic(m1a) # best!

m2 = fit(MixedModel, form2, df, Bernoulli())
m2a = fit(MixedModel, form2a, df, Bernoulli()) # best! but correlation of -1 of random effects!
#m2aa = fit(MixedModel, form2aa, df, Bernoulli()) # no high negative correlation, but worse fit
#m2aaa = fit(MixedModel, form2aaa, df, Bernoulli()) # still high negative correlation

m2b = fit(MixedModel, form2b, df, Bernoulli()) # fixed, performs slightly worse than m2a
m2c = fit(MixedModel, form2c, df, Bernoulli())

bic(m2)
bic(m2a)
bic(m2b)
bic(m2c)

# retro
m3 = fit(MixedModel, form3, df, Bernoulli())
m3a = fit(MixedModel, form3a, df, Bernoulli())
m3b = fit(MixedModel, form3b, df, Bernoulli())

bic(m3) # better than m2a!
bic(m3a)

# scale
m4 = fit(MixedModel, form4, df, Bernoulli()) # not better than m2a!
m4a = fit(MixedModel, form4a, df, Bernoulli()) # better than m2a!
m4b = fit(MixedModel, form4b, df, Bernoulli()) # better than m2a!


bic(m4)
bic(m4a)

# MpD
m5 = fit(MixedModel, form5, df, Bernoulli()) # not better than m4a!
m5a = fit(MixedModel, form5a, df, Bernoulli()) # not better than m4a!

bic(m5)
bic(m5a)

# all
m6 = fit(MixedModel, form6, df, Bernoulli()) # not better than m4a!

