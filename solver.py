#!/usr/bin/python

from random import randrange


def main():
    n = 10
    k = 5

    sol = [randrange(0, n) for _ in range(k)]
    Print(" " * 18 + "Solution:", sol)

    guess = [randrange(0, n) for _ in range(k)]

    info = compare(guess, sol)
    Print(
        f"Possible: {n**k:8d} > Guess:",
        guess,
        f" > c={info[0]} m={info[1]} w={info[2]}",
    )

    possible_seqs = UniversalSet()

    n_steps = 1
    while info[0] != k:
        possible_seqs = generate_seqs(guess, n, info) & possible_seqs
        guess = possible_seqs.pop()
        info = compare(guess, sol)
        n_steps += 1
        Print(
            f"Possible: {1+len(possible_seqs):8d} > Guess:",
            guess,
            f" > c={info[0]} m={info[1]} w={info[2]}",
        )
    Print(" " * 21 + f"Solved in {n_steps} steps!", [])

    return


def compare(g, s):
    k = len(g)

    g_i = [i for i in range(k - 1, -1, -1) if g[i] == s[i]]
    s_j = g_i.copy()

    correct = len(g_i)

    misplaced = 0
    for i in range(k):
        for j in range(k):
            if i not in g_i and j not in s_j and g[i] == s[j]:
                g_i.append(i)
                s_j.append(j)
                misplaced += 1
                break

    wrong = len(g) - correct - misplaced
    return correct, misplaced, wrong


def generate_seqs(g, n, info):
    k = sum(info)

    info_syms = [-3 + i for i, v in enumerate(info) for j in range(v)]
    info_seqs = permutations(info_syms)

    G = set()
    for p in info_seqs:
        G |= sequences(g, p, n)

    return G


def sequences(g, p, n):
    k = len(g)

    M_i = [i for i in range(k) if p[i] == -2]
    M_c = [g[i] for i in M_i]
    W_i = [i for i in range(k) if p[i] == -1]
    W_c = [g[i] for i in W_i]
    avail_colors = [c for c in range(n) if c not in W_c]
    U_i = M_i + W_i

    G = permutations(M_c + [-1] * len(W_c), repel=len(M_c))

    G = complete_sequences(G, U_i, g, p)
    G = expand_placeholders(G, avail_colors, len(W_c))

    return G


def permutations(p0, repel=0):
    p = sorted(p0)

    P = {tuple(p)}

    while 1:
        for j in range(1, len(p)):
            if p[j - 1] < p[j]:
                break
        else:
            return P

        for i in range(j):
            if p[i] < p[j]:
                break

        p[i], p[j] = p[j], p[i]

        p[:j] = reversed(p[:j])

        for r in range(repel):
            if p[r] == p0[r]:
                break
        else:
            P.add(tuple(p))

    return


def complete_sequences(V, inds, g0, p):
    G = set()
    for v in V:
        g = [None] * len(g0)
        for i, j in enumerate(inds):
            g[j] = v[i]
        for i, k in enumerate(p):
            if k == -3:
                g[i] = g0[i]
        G.add(tuple(g))
    return G


def expand_placeholders(G, ac, n_ph):
    while n_ph:
        T = set()
        T, G = G, T
        for t in T:
            t = list(t)
            for i in range(len(t)):
                if t[i] == -1:
                    break
            for j, c in enumerate(ac):
                t[i] = ac[j]
                G.add(tuple(t))
        n_ph -= 1
    return G


def valid(p, g, info):
    new_info = compare(p, g)
    return all([new_info[i] == info[i] for i in range(3)])


class UniversalSet(set):
    def __and__(self, other):
        return other

    def __rand__(self, other):
        return other


def Print(p, s, *args):
    print(p, end=" ")
    print(" ".join([str(t) for t in s]), end="")
    for a in args:
        print(a, end=" ")
    print()
    return


main()
