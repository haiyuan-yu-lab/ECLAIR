import numpy as np
import sys, time, math
from scipy.spatial.distance import pdist,squareform

'''
% Direct Coupling Analysis (DCA)
%
% function dca(inputfile , outputfile)
%
% INPUTS:
%   inputfile  - file containing the FASTA alignment
%   outputfile - file for dca results. The file is composed by N(N-1)/2
%                (N = length of the sequences) rows and 4 columns:
%                residue i (column 1), residue j (column 2),
%                MI(i,j) (Mutual Information between i and j), and
%                DI(i,j) (Direct Information between i and j).
%                Note: all insert columns are removed from the alignment.
%
% SOME RELEVANT VARIABLES:
%   N        number of residues in each sequence (no insert)
%   M        number of sequences in the alignment
%   Meff     effective number of sequences after reweighting
%   q        equal to 21 (20 aminoacids + 1 gap)
%   align    M x N matrix containing the alignmnent
%   Pij_true N x N x q x q matrix containing the reweigthed frequency
%            counts.
%   Pij      N x N x q x q matrix containing the reweighted frequency
%            counts with pseudo counts.
%   C        N(q-1) x N(q-1) matrix containing the covariance matrix.
%
%
% Copyright for this implementation:
%             2011/12 - Andrea Pagnani and Martin Weigt
%                       andrea.pagnani@gmail.com
%                       martin.weigt@upmc.fr
%
% Permission is granted for anyone to copy, use, or modify this
% software and accompanying documents for any uncommercial
% purposes, provided this copyright notice is retained, and note is
% made of any changes that have been made. This software and
% documents are distributed without any warranty, express or
% implied. All use is entirely at the user's own risk.
%
% Aaron Rumack (Haiyuan Yu Lab, Cornell University, 2017) modified the
% DCA algorithm to make it more efficient and translated it into Python.
%
% Any publication resulting from applications of DCA should cite:
%
%     F Morcos, A Pagnani, B Lunt, A Bertolino, DS Marks, C Sander,
%     R Zecchina, JN Onuchic, T Hwa, M Weigt (2011), Direct-coupling
%     analysis of residue co-evolution captures native contacts across
%     many protein families, Proc. Natl. Acad. Sci. 108:E1293-1301.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''

def letter2number(a):
    try:
        return '-ACDEFGHIKLMNPQRSTVWY'.index(a)
    except:
        return 0


def return_alignment(inputfile):
    # reads alignment from inputfile, removes inserts and converts into numbers
    align_full = []
    acc = ''
    for l in open(inputfile):
        if l[0] == '>':
            if acc != '':
                align_full.append(acc)
            acc = ''
        else:
            acc = acc + l.strip()
    align_full.append(acc)
    M = len(align_full)
    ind = map(lambda ch: ch != '.' and ch.upper()==ch,align_full[0])
    N = sum(ind)
    Z = np.zeros((M,N),dtype=int)

    for i in range(M):
        counter = 0
        for j in range(N):
            if ind[j]:
                Z[i,counter] = letter2number(align_full[i][j])
                counter += 1
    q = np.max(Z) + 1
    return [N,M,q,Z]


def get_pi(align,M,N,q,theta):
    W = np.ones((1,M),dtype=np.float32)
    if theta > 0.0:
        W = 1. / (1+np.sum(squareform(pdist(align,'hamming')<theta),axis=0))
    Meff = np.sum(W)

    Pi_true = np.zeros((N,q),dtype=np.float32)

    for l in range(M):
        for k in range(N):
            Pi_true[k,align[l,k]] += W[l]
    Pi_true /= Meff

    return [Pi_true,W]


def pi_pc(Pi_true,pseudocount_weight,N,q):
    return (1.-pseudocount_weight)*Pi_true + pseudocount_weight/q;


def pij_pc(Pij_true, i, j, pseudocount_weight, N, q):
    if i != j:
        return (1.-pseudocount_weight)*Pij_true + (pseudocount_weight/q)/q;
    else:
        scra = np.eye(q);
        return (1.-pseudocount_weight)*Pij_true + (pseudocount_weight/q)*scra;


def true_counts(align,M,q,W,Pi_true,i,j):
    # Returns a q x q matrix that was Pij_true(i,j)
    meff = sum(W)
    Pij_true = np.zeros((q,q),dtype=np.float32)

    if i == j:
        scra = np.eye(q,dtype=np.float32)
        Pij_true = Pi_true[i,:] * scra;
    else:
        for l in range(M):
            Pij_true[align[l,i],align[l,j]] += W[l]
        Pij_true /= meff
    return Pij_true


def mapkey(i,alpha,q):
    return (q-1)*i+alpha


def compute_C(align,M,N,W,pseudocount_weight,Pi_true,Pi,q):
    C = np.zeros((N*(q-1),N*(q-1)),dtype=np.float32)
    for i in range(N):
        for j in range(N):
            pij = pij_pc(true_counts(align,M,q,W,Pi_true,i,j),i,j,pseudocount_weight,N,q)
            alphas = np.repeat(np.arange(q-1),q-1)
            betas = np.tile(np.arange(q-1),q-1)
            idxs1 = mapkey(i,alphas,q)
            idxs2 = mapkey(j,betas,q)
            C[idxs1,idxs2] = pij[alphas,betas] - Pi[i,alphas]*Pi[j,betas]
    print C.dtype
    return C


def calculate_mi(i,j,P2,P1,q):
    M = 0.
    for alpha in range(q):
        for beta in range(q):
            if P2[alpha,beta] > 0:
                M += P2[alpha,beta] * math.log(P2[alpha,beta] / P1[i,alpha]/P1[j,beta])
    return M


def returnW(invC,i,j,q):
    W = np.ones((q,q),dtype=np.float32)
    idxs1 = np.repeat(mapkey(i,np.arange(q-1),q),q-1)
    idxs2 = np.tile(mapkey(j,np.arange(q-1),q),q-1)
    W[:(q-1),:(q-1)] = np.exp(-invC[idxs1,idxs2]).reshape((q-1,q-1))
    return W


def compute_mu(i,j,W,P1,q):
    epsilon = 1e-4
    diff = 1.
    mu1 = np.ones((1,q),dtype=np.float32) / q
    mu2 = np.ones((1,q),dtype=np.float32) / q
    pi = P1[i,:]
    pj = P1[j,:]

    while diff > epsilon:
        scra1 = np.dot(mu2,np.transpose(W))
        scra2 = np.dot(mu1,W)

        new1 = pi / scra1
        new1 = new1/np.sum(new1)

        new2 = pj / scra2
        new2 = new2/np.sum(new2)

        diff = max(np.max(np.abs(new1-mu1)),np.max(np.abs(new2-mu2)))

        mu1 = new1
        mu2 = new2
    return [mu1,mu2]


def compute_di(i,j,W,mu1,mu2,P1):
    tiny = 1.0e-100
    Pdir = W * (np.dot(np.transpose(mu1),mu2))
    Pdir = Pdir / np.sum(Pdir)
    Pfac = np.dot(np.transpose(P1[[i],:]),P1[[j],:])
    return np.trace(np.dot(np.transpose(Pdir),np.log((Pdir+tiny)/(Pfac+tiny))))


def bp_link(i,j,W,P1,q):
    [mu1,mu2] = compute_mu(i,j,W,P1,q)
    return compute_di(i,j,W,mu1,mu2,P1)


def compute_Results(align,M,N,W,Pi,Pi_true,invC,q,fp):
    for i in range(N-1):
        for j in range(i+1,N):
            # Mutual information
            MI_true = calculate_mi(i,j,true_counts(align,M,q,W,Pi_true,i,j),Pi_true,q)

            # Direct information from mean-field
            W_mf = returnW(invC,i,j,q)
            DI_mf_pc = bp_link(i,j,W_mf,Pi,q)
            fp.write('%d %d %g %g\n' % (i+1,j+1,MI_true,DI_mf_pc))
'''
inputfile = sys.argv[1]
outfile = sys.argv[2]

starttime = time.time()
pseudocount_weight = 0.5
theta = 0.2
[N,M,q,align] = return_alignment(inputfile)
[Pi_true,W] = get_pi(align,M,N,q,theta)

print '### N = %d M = %d Meff = %0.2f q = %d\n' % (N,M,sum(W),q)

Pi = pi_pc(Pi_true, pseudocount_weight, N, q)
C = compute_C(align, M, N, W, pseudocount_weight, Pi_true, Pi, q)
print 'Computed C'
invC = np.linalg.inv(C)
print 'Computed C Inverse'

fp = open(outfile,'w')
compute_Results(align,M,N,W,Pi,Pi_true,invC,q,fp)
fp.close()

interval = time.time() - starttime
print 'Time elapsed: %0.f seconds' % interval
'''
