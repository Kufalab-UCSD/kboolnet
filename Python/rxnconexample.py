# Short script with some example stuff on a few of the rxncon Python library's functions
# states and components will be in random orders so just rerun this if you want to see more examples
# NOTE: Some of the terms here are kind confusing, so here's an attempt to explain what each one is
#   - spec: A spec is a "molecule specification", and can describe either an entire component (eg CCR2) or
#           specific domains/residues within that component (eg CCR2_[lig])
#   - component: A component will generally be either a molecule or a protein. Specs can describe components
#   - states: States describe the various states a component can be in. These can be either modifications (eg p38_[(P)]-{p}),
#             interactions with other componenets (eg Gai_[Gbg]--Gbg_[Gai]) or self interactions (idk how those work though)
from rxncon.input.excel_book.excel_book import ExcelBook
from rxncon.visualization.regulatory_graph import RegulatoryGraph
from rxncon.core.state import state_from_str

excel_book = ExcelBook("/home/adrian/Internship/CCR2_rxncon/CCR2_master/out/master.xlsx")
rxncon_system = excel_book.rxncon_system

# The states and components lists hold objects of type Spec
states = rxncon_system.states
components = rxncon_system.components() # For some reason you have to call this method instead of it just being a list

# Shows what a state Spec and a component Spec look like
print(states[0])
print(components[0])

# You can pick these out by name using their "name" attribute (list comprehension woo)
print([state for state in states if (state.name == "CCR2_[Ga]--0" or state.name == "Pyk2_[(P)]-{{p}}")])
print([component for component in components if (component.name == "AC" or component.name == "CCL2")])

# You can see which states "live" on a component (i.e. which states involve a given component)
print(rxncon_system.states_for_component(components[0]))
print(rxncon_system.states_for_component_grouped(components[0])) # This version outputs a dict grouped by mutually exclusive states

# You can see which components a state involves
print(states[0].components)

# You can also see the mutually exclusive states for a given state
print(rxncon_system.complement_states(states[0]))

# You can get fully expanded state(s) from a string
orig_state = state_from_str("Gai--0") 
expanded_states = [state for state in states if state.is_subset_of(orig_state)]
print(expanded_states)

# You can get the reactions that produce a given state
producing_reactions = [reaction for reaction in rxncon_system.reactions if any([state.is_subset_of(orig_state) for state in reaction.produced_states])]
print(producing_reactions)

# You can see which reactions are bidirectional using the BIDIRECTIONAL_REACTIONS variable
from rxncon.core.reaction import BIDIRECTIONAL_REACTIONS
print(BIDIRECTIONAL_REACTIONS)

# There's definitely more stuff in the rxncon source files, let me know if you want me to try and dig some new stuff up!
