using SpectralFitting
using Test

include("../dummies.jl")

# single model
m = PowerLaw()
m2 = copy(m)
@inferred copy(m)
m.K.value = 2.0
@test m2.K.value == 1.0

# table model
m = DummyMultiplicativeTableModel()
m2 = copy(m)
@inferred copy(m)
m.a.value = 2.0
@test m2.a.value == 1.0
@test m.table === m2.table

# composite model
model = PowerLaw() + PowerLaw()
model2 = copy(model)
model.K_1.value = 2.0
@test model2.K_2.value == 1.0

@inferred copy(model)
