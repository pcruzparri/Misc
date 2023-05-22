import WrightTools as wt
import attune
from attune._transition import Transition
from attune._tune import Tune

"""
SCRIPT NOTES HERE

"""


instrument_name = 'w3'
instrument = attune.load(instrument_name)####### 
arr_from = 'sig'
arr_to = 'idl'
tune = 'bbo' 
degeneracy = 1600 #in nm 
degen_units = 'nm'
instr_port = 39303
apply = False
md = {
    "arrangement": arr_to,
    "tune": tune,
    "degeneracy": degeneracy,
    "degen_units": degen_units,
}        

# Get all the previous arrangements in the instrument
arrs = {} 
for arr in instrument.as_dict()['arrangements'].keys():
    arrs[arr] = instrument.as_dict()['arrangements'][arr]

# Generate the new independent points from prev arrangement to new arrangement    
from_tune = instrument[arr_from][tune] 

if arr_from == 'idl' and arr_to = 'sig':
    idl_points = from_tune.independent
    new_independents = [idl + 2*(amount-idl) for idl in idl_points]
if arr_from == 'sig' and arr_to = 'idl':
    sig_points = from_tune.independent
    new_independents = [sig - 2*(sig-degeneracy) for sig in sig_points]

arrs[arr_to] = attune.Arrangement(arr_to, dict(tune=from_tune))     

# Generate new instrument with new arrangement included 
instr = attune.Instrument(arrs, name=instrument_name)

if amount_units is not None:
    amount = wt.units.convert(amount, amount_units, from_tune.dep_units)

# Generate the new tune for the new arrangement
instr[arr_to]._tunes[tune] = Tune(
    new_independents,
    from_tune.dependent,
    dep_units=from_tune.dep_units,
)

# Transition to new instrument with new tuned arrangement
instr._transition = Transition(f'adding_{arr_from}_to_{arr_to}', instrument, metadata=md)
instr._load = None

# Apply the new instrument to the client
if apply:
    import yaqc
    c = yaqc.Client(instr_port)
    c.set_instrument(instr.as_dict())

#call the new instrument object for proof checking on terminal
instr


