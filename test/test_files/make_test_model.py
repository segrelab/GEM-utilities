from cobra import Metabolite, Model, Reaction

model = Model("example_model")

reaction = Reaction("R_3OAS140")
reaction.name = "3 oxoacyl acyl carrier protein synthase n C140 "
reaction.subsystem = "Cell Envelope Biosynthesis"
reaction.lower_bound = 0.0  # This is the default
reaction.upper_bound = 1000.0  # This is the default

ACP_c = Metabolite(
    "ACP_c", formula="C11H21N2O7PRS", name="acyl-carrier-protein", compartment="c"
)
omrsACP_c = Metabolite(
    "M3omrsACP_c",
    formula="C25H45N2O9PRS",
    name="3-Oxotetradecanoyl-acyl-carrier-protein",
    compartment="c",
)
co2_c = Metabolite("co2_c", formula="CO2", name="CO2", compartment="c")
malACP_c = Metabolite(
    "malACP_c",
    formula="C14H22N2O10PRS",
    name="Malonyl-acyl-carrier-protein",
    compartment="c",
)
h_c = Metabolite("h_c", formula="H", name="H", compartment="c")
ddcaACP_c = Metabolite(
    "ddcaACP_c",
    formula="C23H43N2O8PRS",
    name="Dodecanoyl-ACP-n-C120ACP",
    compartment="c",
)

reaction.add_metabolites(
    {malACP_c: -1.0, h_c: -1.0, ddcaACP_c: -1.0, co2_c: 1.0, ACP_c: 1.0, omrsACP_c: 1.0}
)


reaction.gene_reaction_rule = "( STM2378 or STM1197 )"

model.add_reactions([reaction])

model.objective = "R_3OAS140"
