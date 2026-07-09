
import random
from collections import defaultdict, deque
from typing import List, Tuple, Dict, Set

Individual = Tuple[List[int], List[int]]  # (seq, x)

# ---------------------------
# Utilities for precedence
# ---------------------------

def topo_sort(seq: List[int], precedence: Dict[int, List[int]]) -> List[int]:
    """
    Repair helper: returns a topologically valid permutation using
    the relative order in 'seq' as a tie-breaker to preserve intent.
    """
    # Build graph / indegree from precedence restricted to items in seq
    present = set(seq)
    preds = {j: [p for p in precedence.get(j, []) if p in present] for j in present}
    indeg = {j: 0 for j in present}
    for j, ps in preds.items():
        for p in ps:
            indeg[j] += 1

    # stable queue using original order for determinism
    order_index = {op: i for i, op in enumerate(seq)}
    q = [j for j in present if indeg[j] == 0]
    q.sort(key=lambda j: order_index[j])

    out = []
    while q:
        j = q.pop(0)
        out.append(j)
        for k in present:
            if j in preds.get(k, []):
                indeg[k] -= 1
                if indeg[k] == 0:
                    q.append(k)
                    q.sort(key=lambda x: order_index[x])
    if len(out) != len(seq):
        # Cycle or missing items detected; as a fallback keep original order
        return seq[:]
    return out


def ensure_flags_feasible(seq: List[int], x: List[int],
                          precedence: Dict[int, List[int]]) -> List[int]:
    """
    Ensures that if an operation is marked performed (x[i]==1),
    all its predecessors in 'seq' are performed as well when required
    by your domain logic. Here we adopt a conservative rule:
    a performed operation implies all listed predecessors must be performed.
    """
    pos = {op: i for i, op in enumerate(seq)}
    x = x[:]
    for i, op in enumerate(seq):
        if x[i] == 1:
            for p in precedence.get(op, []):
                if p in pos:
                    x[pos[p]] = 1  # enforce predecessor execution
    return x


# ---------------------------
# Representation helpers
# ---------------------------

def make_individual(seq: List[int], x: List[int],
                    precedence: Dict[int, List[int]]) -> Individual:
    """
    Create/repair an individual to be feasible:
    - topologically sort seq under precedence
    - adjust flags x to be consistent with precedence
    """
    assert len(seq) == len(x), "seq and x must be the same length"
    # 1) Remove duplicates while preserving original order (parents may overlap)
    seen: Set[int] = set()
    seq_unique, x_unique = [], []
    for op, flag in zip(seq, x):
        if op not in seen:
            seen.add(op)
            seq_unique.append(op)
            x_unique.append(1 if flag else 0)
    # 2) Topological repair
    seq_fixed = topo_sort(seq_unique, precedence)
    # 3) Re-align flags with new order
    index_from_old = {op: i for i, op in enumerate(seq_unique)}
    x_fixed = [x_unique[index_from_old[op]] for op in seq_fixed]
    # 4) Enforce feasibility of flags
    x_fixed = ensure_flags_feasible(seq_fixed, x_fixed, precedence)
    return seq_fixed, x_fixed


# ---------------------------
# PPX crossover (paper Fig. 2)
# ---------------------------

def ppx_crossover(parent1: Individual, parent2: Individual,
                  precedence: Dict[int, List[int]],
                  rng: random.Random = random) -> Individual:
    """
    Precedence-Preserving Crossover (PPX)
    Matches the paper's pseudocode:
      - At each step pick next available item from P1 or P2 with prob 0.5
      - Copy the corresponding execution flag from the same parent
      - Remove that front item from that parent
      - Continue until exhausted
    Finally, repair feasibility against precedence.
    """
    p1_seq, p1_x = list(parent1[0]), list(parent1[1])
    p2_seq, p2_x = list(parent2[0]), list(parent2[1])

    child_seq, child_x = [], []

    # For convenience, we consume from the "front". Use lists as queues.
    while p1_seq or p2_seq:
        choose_p1 = False
        if not p2_seq:
            choose_p1 = True
        elif not p1_seq:
            choose_p1 = False
        else:
            choose_p1 = (rng.random() < 0.5)

        if choose_p1:
            op = p1_seq.pop(0)
            flag = p1_x.pop(0)
        else:
            op = p2_seq.pop(0)
            flag = p2_x.pop(0)

        # If operation already added, skip (needed because parents share ops)
        if op in child_seq:
            continue

        child_seq.append(op)
        child_x.append(1 if flag else 0)

        # Also remove this op from the front of the other parent if it happens to be next
        # (not required, but keeps behavior close to the pseudocode's "delete" action)
        if p1_seq and p1_seq[0] == op:
            p1_seq.pop(0); p1_x.pop(0)
        if p2_seq and p2_seq[0] == op:
            p2_seq.pop(0); p2_x.pop(0)

    # Feasibility repair as per paper
    return make_individual(child_seq, child_x, precedence)


# ---------------------------
# Mutation operators (paper)
# ---------------------------

def mutate_swap(ind: Individual, precedence: Dict[int, List[int]],
                rng: random.Random = random) -> Individual:
    """Swap two positions in seq, then repair."""
    seq, x = list(ind[0]), list(ind[1])
    if len(seq) >= 2:
        i, j = rng.sample(range(len(seq)), 2)
        seq[i], seq[j] = seq[j], seq[i]
        x[i], x[j] = x[j], x[i]  # keep flags attached to operations
    return make_individual(seq, x, precedence)


def mutate_toggle(ind: Individual, precedence: Dict[int, List[int]],
                  p_flip: float = 0.15, rng: random.Random = random) -> Individual:
    """Randomly toggle execution flags, then repair."""
    seq, x = list(ind[0]), list(ind[1])
    for i in range(len(x)):
        if rng.random() < p_flip:
            x[i] = 1 - x[i]
    return make_individual(seq, x, precedence)


def mutate(ind: Individual, precedence: Dict[int, List[int]],
           rng: random.Random = random) -> Individual:
    """Randomly choose one mutation operator, as in the paper."""
    if rng.random() < 0.5:
        return mutate_swap(ind, precedence, rng)
    else:
        return mutate_toggle(ind, precedence, rng)


# ---------------------------
# Monte-Carlo objective evaluation (stub)
# ---------------------------

def evaluate_objectives_monte_carlo(
    ind: Individual,
    samples: int,
    draw_times,      # callable: op_id -> float (disassembly time sample)
    draw_setup,      # callable: (prev_op, op) -> float (setup time sample)
    profit_of,       # callable: op_id -> float (recycling/reuse value)
    cost_rate_of,    # callable: op_id -> float (cost per time unit)
    setup_cost_rate, # callable: (prev_op, op) -> float (cost per time unit in setup)
    energy_rate_of,  # callable: op_id -> float (energy per time unit)
    setup_energy_rate # callable: (prev_op, op) -> float (energy per time unit in setup)
):
    """
    Returns (E[profit], E[energy]) estimated by Monte-Carlo as in the paper.
    You can plug this into NSGA-II or MOEA/D.
    """
    seq, x = ind
    # Build the executed subsequence in order
    executed = [op for op, flag in zip(seq, x) if flag == 1]

    total_profit, total_energy = 0.0, 0.0
    for _ in range(samples):
        profit = 0.0
        energy = 0.0
        prev = None
        for op in executed:
            # profit from subassembly (reuse/recycle value)
            profit += profit_of(op)

            # disassembly operation time -> cost & energy
            t = draw_times(op)
            c_rate = cost_rate_of(op)
            e_rate = energy_rate_of(op)
            # profit net of cost (paper has profit minus cost terms)
            profit -= c_rate * t
            energy += e_rate * t

            # setup if applicable
            if prev is not None:
                ts = draw_setup(prev, op)
                cs_rate = setup_cost_rate(prev, op)
                es_rate = setup_energy_rate(prev, op)
                profit -= cs_rate * ts
                energy += es_rate * ts
            prev = op

        total_profit += profit
        total_energy += energy

    # Expected values (the paper uses expectations in the objectives)
    return total_profit / samples, total_energy / samples
