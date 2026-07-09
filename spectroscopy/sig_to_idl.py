import copy

from attune._transition import Transition
from attune._tune import Tune


def sig_to_idl(instrument, arrangement, tune, amount, amount_units=None): # MUST PASS SIG AS TUNE AND DEGENERACY AS AMOUNT. #amount in nm  
    md = {
        "arrangement": arrangement,
        "tune": tune,
        "amount": amount,
        "amount_units": amount_units,
    }
    sig_tune = instrument[arrangement][tune]
    idler = attune.Arrangement("idler", dict(bbo=sig_tune))
    def _get_sig_offset(sig_tune, amount):
        sigs = sig_tune.independent
        return [sig - 2*(sig-amount) for sig in sigs]
    if amount_units is not None:
        amount = wt.units.convert(amount, amount_units, sig_tune.dep_units)
    instr = copy.deepcopy(instrument)
    instr["idler"]._tunes[tune] = Tune(
        _get_sig_offset,
        sig_tune.dependent,
        dep_units=sig_tune.dep_units,
    )
    instr._transition = Transition("offset_by", instrument, metadata=md)
    instr._load = None
    return instr