import numpy as np
from itertools import groupby

class dpmeans:

    def __init__(self, X, w=0, cl_graph=None, Lambda=None, scores=0):
        # Initialize parameters for DP means
        self.K = 1
        self.K_init = 4
        self.d = X.shape[1]
        self.z = np.mod(np.random.permutation(X.shape[0]),self.K)+1
        self.mu = np.random.standard_normal((self.K, self.d))
        self.sigma = 1
        self.nk = np.zeros(self.K)
        self.pik = np.ones(self.K)/self.K 
        
        self.cl_graph = cl_graph # cannot link elements
        self.w = w # penalty of cl elements
        self.scores = np.array(scores) # binding scores
        
        #init mu
        self.mu = np.array([np.mean(X,0)])
        
        #init lambda
        if not Lambda:
            self.Lambda = self.kpp_init(X,self.K_init)
        else:
            self.Lambda = Lambda
        
        self.max_iter = 500
        self.obj = np.zeros(self.max_iter)
        self.em_time = np.zeros(self.max_iter)   
        
    def kpp_init(self,X,k):
        #k++ init
        #lambda is max distance to k++ means

        [n,d] = np.shape(X)
        mu = np.zeros((k,d))        
        dist = np.inf*np.ones(n)
                
        mu[0,:] = X[int(np.ceil(np.random.rand()*n-1)),:]
        for i in range(1,k):
            D = X-np.tile(mu[i-1,:],(n,1))
            dist = np.minimum(dist, np.sum(D*D,1))
            idx = np.where(np.random.rand() < np.cumsum(dist/float(sum(dist))))
            mu[i,:] = X[idx[0][0],:]
            Lambda = np.max(dist)
        
        return Lambda

    def _objective_function(self, X, x_i, dist, z):
        d = dist[x_i]
        
        penalty = np.zeros(self.K)
        if self.cl_graph:
            cl_z = np.array([z[y_i] for y_i in self.cl_graph[x_i]])
            for i in cl_z[cl_z!=-1]:
                penalty[i] += self.w
        return d+penalty
    
    def fit(self,X):

        obj_tol = 1
        max_iter = self.max_iter        
        [n,d] = np.shape(X)
        
        obj = np.zeros(max_iter)
        em_time = np.zeros(max_iter)
#         print('running dpmeans...')
        
        for iter in range(max_iter):
            dist = np.zeros((n,self.K))
            
            #assignment step
            for kk in range(self.K):
                Xm = X - np.tile(self.mu[kk,:],(n,1))
                dist[:,kk] = np.sum(Xm*Xm,1)
            
            #update labels
            dmin = np.min(dist,1)
            
            z = np.full(n, fill_value=-1)
            index = list(range(X.shape[0]))
            np.random.shuffle(index)
            for x_i in index:
                z[x_i] = np.argmin(self._objective_function(X, x_i, dist, z))
                
                if dmin[x_i] > self.Lambda:
                    self.K = self.K + 1
                    z[x_i] = self.K-1 #cluster labels in [0,...,K-1]
                    self.mu = np.vstack([self.mu, X[x_i].reshape([1,-1])])
                    Xm = X - np.tile(self.mu[self.K-1,:],(n,1))
                    dist = np.hstack([dist, np.array([np.sum(Xm*Xm,1)]).T])
                    dmin = np.min(dist,1)
                
            self.z = z
            
            #update step
            self.nk = np.zeros(self.K)
            delete = []
            for kk in range(self.K):
                self.nk[kk] = self.z.tolist().count(kk)
                if self.nk[kk]==0:
                    self.K -= 1
                    delete.append(kk)
                else:
                    idx = np.where(self.z == kk)
                    self.mu[kk,:] = np.mean(X[idx[0],:],0)
                    if self.scores.any() != 0:
                        start = min(max(0, int(self.mu[kk,:]-35)), len(self.scores))
                        end = min(max(0, int(self.mu[kk,:]+35)), len(self.scores))
                        self.mu[kk,:] = start+np.argmin((abs(X[idx[0],:]-np.arange(start,end)-self.scores[start:end])).sum(0))

            self.mu = np.delete(self.mu, delete, 0)
            
            self.pik = self.nk/float(np.sum(self.nk))
            
            #compute objective
            for kk in range(self.K):
                idx = np.where(self.z == kk)
                obj[iter] = obj[iter] + np.sum(dist[idx[0],kk],0)          
            obj[iter] = obj[iter] + self.Lambda * self.K
            
            #check convergence
            if (iter > 0 and np.abs(obj[iter]-obj[iter-1]) < obj_tol*obj[iter]):
                break
        return self.z, self.mu

def encode_seq(seq):
    seq = np.array([seq[i:i+2] for i in range(0, len(seq), 2)])
    return np.where((seq == 'AT') | (seq == 'TA') | (seq == 'AA') | (seq == 'TT'), 1, 0)

def fasta_iter(fasta_name):
    """
    modified from Brent Pedersen
    Correct Way To Parse A Fasta File In Python
    given a fasta file. yield tuples of header, sequence
    """
    "first open the file outside "
    fh = open(fasta_name)

    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

    for header in faiter:
        # drop the ">"
        headerStr = header.__next__()[1:].strip()

        # join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())

        yield (headerStr, seq)
