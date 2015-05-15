import sys
sys.ps1 = 'SOMETHING'
from matplotlib import pyplot as plt
import numpy as np
import random
from itertools import cycle
from collections import Counter
import copy
from matplotlib import interactive

# Please run %matplotlib qt before %run StatConnProject.py

class CA3_Model():
    '''EN.580.694 Statistical Connectomics 
    SBM class, including analysis methods
    William Hockeimer, JHU Neuroscience
    Parameters    nverts -- number of vertices
                   kblocks -- number of blocks (assume symmetric graph)
                   baseps -- high and low prob for on/off diag base prob
    '''

    def __init__(self,nverts):
        if nverts%9 != 0:
            return 'Number of vertices must be divisible by 9'
        self.nverts = nverts
        self.kblocks = 9 # 9 expression zones. I know magic numbers are bad. Sue me. 
        self.baseps = (0,0)
        self.obsv_ca3 = {}
        self.obsv_ca3_index = 0
        self.cues = {
            'cdh10': (0,0,1,1,1,0,0,0,0),
            'cdh11': (0,1,1,1,1,0,0,0,0),
            'cdh2': (0,0,1,0,1,0,0,0,0),
            'cdh24': (0,1,1,1,1,0,0,0,0),
            'cdh8': (0,1,1,0,1,0,0,0,0),
            'col12a1': (0,0,0,1,0,0,0,1,0),
            'col15a1': (0,0,1,1,0,0,0,0,0),
            'col6a1': (0,0,0,0,0,0,0,0,1),
            'col5a1': (0,0,0,1,0,0,1,1,1),
            'col6a2': (0,0,0,0,0,0,0,1,1),
            'efna5': (0,0,0,0,0,0,0,1,1),
            'efnb2': (0,0,0,0,0,0,0,1,1),
            'efnb3': (0,0,0,0,0,1,0,1,1),
            'slit2': (0,1,1,1,1,0,0,0,1)
                    }
        self.receptors = {
            'cdh10': (0,0,1,1,1,0,0,0,0),
            'cdh11': (0,1,1,1,1,0,0,0,0),
            'cdh2': (0,0,1,0,1,0,0,0,0),
            'cdh24': (0,1,1,1,1,0,0,0,0),
            'cdh8': (0,1,1,0,1,0,0,0,0),
            'col12a1': (0,0,0,1,0,0,0,1,0),
            'col15a1': (0,0,1,1,0,0,0,0,0),
            'col6a1': (0,0,0,0,0,0,0,0,1),
            'col5a1': (0,0,0,1,0,0,1,1,1),
            'col6a2': (0,0,0,0,0,0,0,1,1),
            'epha7': (1,1,1,1,0,0,0,0,0),
            'ephb1': (0,1,1,0,0,0,0,0,0),
            'epha3': (0,1,1,0,0,0,0,0,0),
            'epha4': (1,1,1,1,0,1,1,0,0),
            'robo1': (1,1,1,1,1,0,0,0,0),
            'robo2': (1,1,1,1,1,0,0,0,0)
                    }
        self.fam_dict = {
            'cad': ['cdh10','cdh11','cdh24','cdh8','cdh2'],
            'coll': ['col12a1','col15a1','col6a1','col6a2','col5a1'],
            'ephA': ['epha7','epha3','epha4','efna5'],
            'ephB': ['efnb2','efnb3','efnb1'],
            'robo': ['slit2','robo1','robo2']
                    }
        
    def draw_sbm(self,cond='ctrl',baseps=(0.75,0.25),maxmod=0.3,plot=True):
        '''Draw an observed CA3 SBM model from a distribution governed by betamat
        PARMS:      cond = ctrl for classic recurrent net, exp for genetically rewired net
                    base = if exp, this sets baseline connection prob
                    maxmod = if exp, this sets maximum amt of rewiring 
                    plot = plot spy of adj mat?
                    
        '''
        #  Define vertex labels and base probs
        self.baseps = baseps
       


        ## Control: Define the beta matrix of inter- and intra-block connection probs
        if cond == 'ctrl':
            betamat = self.gen_ctrl_betamat()
        elif cond == 'exp':
            self.maxmod = maxmod
            betamat = self.gen_betamat()


        ## Define adjacency matrix
        adj,labels = self.define_adj(betamat) 
       
        ## Save adjacency matrix, beta matrix, and node labels
        self.obsv_ca3[self.obsv_ca3_index] = (cond,labels,baseps,betamat,adj)
        self.obsv_ca3_index += 1
        
        ##Plot
        if plot:
            plt.spy(adj)
            if cond == 'exp':
                plt.title('Adj Matrix {}: {} CA3\n base prob: {}, max mod: {}'.format(self.obsv_ca3_index,cond,self.baseps[0],maxmod))
            else:
                plt.title('Adj Matrix {}: {} CA3\n probs: {},{}'.format(self.obsv_ca3_index,cond,self.baseps[0],self.baseps[1]))
            plt.xlabel('From Vertex...',fontsize=16)
            plt.ylabel('To Vertex...',fontsize=16)
            plt.show()

    
    def gen_betamat(self):
        ''' Generate beta matrix using C/R interaction table'''

        betamat_delta = [[[] for i in range(9)] for j in range(9)]
        betamat = [[self.baseps[0] for i in range(9)] for j in range(9)]
        for from_zone in range(self.kblocks):
            for to_zone in range(self.kblocks):
                self.eval_genes(from_zone,to_zone,betamat_delta)
        betamat = self.modulate_betamat(betamat,betamat_delta)
        return betamat
    
    def eval_genes(self,i,j,betamat_delta):
        '''Given pairing of zone i connecting to zone j, evaluate all genetic
        interactions to determine degree of connectivity modulation from baseline'''
        for k,cue in enumerate(self.cues):
            if self.cues[cue][i] == 1:
                fam,famlist = self.get_fam(cue)
                recs_in_fam = [rec for rec in self.receptors if rec in famlist]
                for rec in recs_in_fam:
                    if self.receptors[rec][j] == 1:
                        betamat_delta = self.apply_rule(cue,rec,i,j,fam,betamat_delta)
                        
    def get_fam(self,cue):
        ''' Find family of cue, return string and list of it'''
        for k,v in enumerate(self.fam_dict):
            if cue in self.fam_dict[v]:
                return v,self.fam_dict[v]
                        
                        
    def apply_rule(self,cue,rec,i,j,fam,betamat_delta):
        '''Given a C/R family, apply the relevant rule to determine if there's a hit'''
        if fam == 'coll':
            if cue == rec:
                betamat_delta[i][j].append(1.0)

        if fam == 'cad':
            if cue == rec:
                betamat_delta[i][j].append(1.0)

        if fam == 'ephA':
            betamat_delta[i][j].append(-1.0)

        if fam == 'ephB':
            betamat_delta[i][j].append(-1.0)

        if fam == 'robo':
            betamat_delta[i][j].append(-1.0)
    
        return betamat_delta
    
    def modulate_betamat(self,betamat,betamat_delta,maxmod=0.3):
        ''' Notes'''
        tot = max(self.cues,self.receptors)
        for i in range(self.kblocks):
            for j in range(self.kblocks):
                if len(betamat_delta[i][j]) != 0:
                    betamat[i][j] = betamat[i][j] + (sum(betamat_delta[i][j])/len(betamat_delta[i][j]))*self.maxmod
        return betamat
    
    def define_adj(self,betamat):
        '''Define adj mat using betamat. Constrain number of edges.'''
        
        cpv = (72000*self.nverts)/300000 # Empirically ballparked num. conns/cell, total cells
        vpb = self.nverts/self.kblocks
        labels = [i*np.ones(vpb) for i in range(self.kblocks)]
        labels = [item for sublist in labels for item in sublist]
        labels = [int(i) for i in labels]
        adj = np.asarray([[0 for i in range(1,self.nverts+1)] for j in range(1,self.nverts+1)])
        for v in range(1,self.nverts+1):
            verts_tested = []
            c = cpv
            zcyc = cycle(range(self.kblocks))
            z = zcyc.next()
            while c > 0:
                t = random.randint(1+(z*vpb),vpb+(z*vpb))
                if t not in verts_tested:
                    verts_tested.append(t)
                    if random.random() < betamat[labels[t-1]][labels[v-1]]:
                        adj[t-1][v-1] = 1
                        c += -1
                    z = zcyc.next()
        return adj,labels
    
    def avg_degree(self,sbmidx):
        '''Compute average incoming and outgoing degree for an
        PARMS:    sbmidx - 0-indexed key of sbm to use
        RETUTRNS: avg incoming degee, avg outgoing degree
        '''
        adj = self.obsv_sbm[sbmidx][0]
        in_deg = [float(sum(adj[i] == True)) for i in range(self.nverts)]
        avg_in_deg = np.mean(in_deg)
        out_deg = [float(sum(adj[:,i] == True)) for i in range(self.nverts)]
        avg_out_deg = np.mean(out_deg)
        return avg_in_deg,avg_out_deg
    
    def find_hubs(self,adj,labels):
        '''
        Given a threshold, find nodes which are hubs from
        perspective of incoming and outgoing nodes.
        PARMS:    thresh - num of std devs of edges above which a node is an in/out hub
        '''
        thresh = 1.5
        in_deg = [float(sum(adj[i] == True)) for i in range(self.nverts)]
        out_deg = [float(sum(adj[:,i] == True)) for i in range(self.nverts)]
        thresh_in = (thresh*np.std(in_deg)) + np.mean(in_deg)
        thresh_out = (thresh*np.std(out_deg)) + np.mean(out_deg)
        in_hubs = [i >= thresh_in for i in in_deg]
        out_hubs = [i >= thresh_out for i in out_deg]
        in_hubs = list(in_hubs)
        if len(in_hubs)==0:
            return (0,0)
        hubs_idx = [k for k,v in enumerate(in_hubs) if v == True]
        hubs_labels = [labels[i] for i in hubs_idx]
        c = Counter(hubs_labels)
        c = c.items()
        num_hubs = sum([c[i][1] for i in range(len(c)) ])
        return num_hubs
    
    
    def plot_gene_info(self,fig):
        '''Plot distribution of gene expression across zones'''
        
        cue_mat = [ca3.cues[cue] for k,cue in enumerate(ca3.cues)]
        rec_mat = [ca3.receptors[rec] for k,rec in enumerate(ca3.receptors)]
        cues = [v for k,v in enumerate(self.cues)]
        receptors = [v for k,v in enumerate(self.receptors)]
        zones = ['Zone {}'.format(i+1) for i in range(9)]
        fig.add_subplot(1,2,1)
        plt.spy(np.asarray(cue_mat).T)
        plt.title('Distribution of Axon Guidance Cues in Hippocampal CA3',fontsize=26,y=1.4)
        plt.xticks(range(len(cues)),cues,size='large',rotation='vertical')
        plt.yticks(range(len(zones)),zones,size='large')
        fig.add_subplot(1,2,2)
        plt.spy(np.asarray(rec_mat).T)
        plt.title('Distribution of Axon Guidance Cue Receptors in Hippocampal CA3',fontsize=26,y=1.4)
        plt.xticks(range(len(receptors)),receptors,size='large',rotation='vertical')
        plt.yticks(range(len(zones)),zones,size='large')
        plt.show()
        
    def hub_analysis(self,numiters=100):
        '''Generate many sbms just saving hub info. Analyze across sbm draws.
        Perform permutation test. Similar to self.find_hubs except generating 
        sbms and not saving them, save their hub info.'''
        
        ctrl_adj= self.obsv_ca3[0][4]
        exp_adj = self.obsv_ca3[1][4]
        labels = self.obsv_ca3[1][1]
        ctrl_count = self.find_hubs(ctrl_adj,labels)
        exp_count = self.find_hubs(exp_adj,labels)
        shuff_list = []

        for i in range(numiters):
            shuff_adj = copy.deepcopy(exp_adj)
            self.shuffle_mat(shuff_adj)
            shuff_list.append(self.find_hubs(shuff_adj,labels))

      
        shuff_count = sum(shuff_list)/numiters #normalize
        shuff_hist = np.histogram(shuff_list)
        perm_p  = float(len([i for i in shuff_list if i > exp_count]))/float(len(shuff_list))

        return ctrl_count,exp_count,shuff_count,perm_p,shuff_list

    
    def gen_ctrl_betamat(self):
        ''' Generate am SBM prob distribution corresponding
        to a fully recurrent net (with 2 probs: on and off diag)'''
        betamat_ctrl = np.asarray([[0.0 for i in range(self.kblocks)] for j in range(self.kblocks)])
        for ii in range(self.kblocks):
            for jj in range(self.kblocks):
                if ii == jj:
                    betamat_ctrl[ii][jj] = self.baseps[0]  + 0.15*random.random()
                else:
                    betamat_ctrl[ii][jj] = self.baseps[1]  + 0.15*random.random()
        return betamat_ctrl
    
    def shuffle_mat(self,mat):
        '''Given an adj mat, shuffle it for a permutation test'''
        [np.random.shuffle(mat[i]) for i in range(len(mat))]
        [np.random.shuffle(mat[:,i]) for i in range(len(mat))]
        return mat
        
        
        

        
        
# Please run %matplotlib qt before %run StatConnProject.py
  

        
def routine(model,baseps,maxmod,shuff_counts,mean_exp):
	'''Given an instance of CA3_Model and input parameters below, run
	and graph model'''
	model.draw_sbm('ctrl',baseps,maxmod,False)
	model.draw_sbm('exp',baseps,maxmod,False)
	fig1 = plt.figure()
	fig1.add_subplot(2,2,1)
	plt.spy(model.obsv_ca3[0][4])
	plt.ylabel('To Vertex...',fontsize=20)
	plt.xlabel('From Vertex...',fontsize=20)
	plt.title('Control CA3 Adjacency Matrix',fontsize=24)
	fig1.add_subplot(2,2,2)
	plt.spy(model.obsv_ca3[1][4])
	plt.ylabel('To Vertex...',fontsize=18)
	plt.xlabel('From Vertex...',fontsize=18)
	plt.title('Genetically Re-wired CA3 Adjacency Matrix',fontsize=24)
	fig1.add_subplot(2,2,3)
	plt.pcolor(model.obsv_ca3[0][3])
	plt.colorbar()
	plt.title('Control Graph Beta Matrix',fontsize=24)
	fig1.add_subplot(2,2,4)
	plt.pcolor(np.asarray(model.obsv_ca3[1][3]))
	plt.colorbar()
	plt.title('Experimental Graph Beta Matrix',fontsize=24)
	
	fig2 = plt.figure()
	model.plot_gene_info(fig2)
	
	fig3 = plt.figure()
	plt.hist(shuff_counts)
	plt.vlines(mean_exp,0,80)
	plt.title('Observed Hubs vs. Hubs from 500 Shuffled Graphs',fontsize=26)
	plt.xlabel('Number of Hubs',fontsize=20)
	plt.ylabel('Frequency of Occurance',fontsize=20)
	plt.show()
	
shuff_counts = [61, 68, 59, 69, 62, 67, 60, 60, 62, 69, 66, 77, 61, 71, 69, 64, 72, 64, 63, 70, 76, 71, 62, 75, 68, 63, 77, 75, 80, 71, 77, 60, 74, 64, 77, 70, 66, 71, 60, 87, 73, 70, 69, 66, 64, 76, 67, 78, 75, 66, 79, 72, 
63, 63, 64, 68, 72, 65, 72, 68, 62, 65, 77, 72, 55, 73, 66, 62, 59, 67, 60, 66, 56, 65, 77, 56, 70, 64, 63, 68, 60, 65, 82, 73, 64, 71, 62, 64, 69, 59, 65, 73, 72, 70, 71, 62, 72, 67, 75, 70, 61, 65, 78, 65, 70, 66, 57, 82, 
67, 70, 73, 68, 78, 72, 68, 70, 68, 71, 63, 72, 69, 72, 72, 67, 52, 69, 61, 81, 66, 73, 73, 66, 68, 63, 77, 63, 73, 75, 60, 64, 64, 72, 68, 60, 57, 74, 58, 79, 66, 64, 67, 66, 71, 67, 66, 70, 71, 68, 61, 68, 73, 62, 74, 74, 
76, 74, 57, 77, 66, 73, 63, 71, 63, 71, 71, 75, 63, 61, 69, 66, 72, 79, 59, 72, 67, 66, 76, 73, 55, 55, 58, 57, 75, 63, 74, 68, 69, 69, 72, 78, 75, 68, 82, 66, 66, 70, 65, 60, 82, 84, 62, 64, 72, 66, 64, 62, 75, 68, 77, 74, 
74, 62, 79, 74, 70, 67, 61, 74, 67, 67, 74, 59, 59, 68, 65, 69, 71, 73, 68, 77, 72, 71, 63, 73, 62, 64, 63, 68, 66, 64, 72, 68, 66, 75, 65, 71, 72, 62, 65, 72, 70, 73, 70, 63, 67, 60, 77, 78, 73, 67, 72, 76, 60, 71, 74, 63, 
78, 81, 60, 68, 67, 74, 66, 63, 66, 66, 74, 66, 67, 70, 78, 58, 62, 61, 70, 71, 61, 61, 67, 60, 79, 66, 64, 71, 64, 60, 63, 76, 71, 81, 70, 64, 69, 60, 62, 62, 65, 71, 73, 73, 64, 67, 78, 68, 65, 58, 76, 73, 70, 75, 70, 61, 
63, 64, 83, 67, 57, 67, 70, 72, 54, 70, 69, 70, 71, 81, 69, 60, 74, 66, 65, 73, 66, 73, 71, 76, 69, 69, 70, 80, 73, 68, 70, 76, 59, 75, 68, 64, 60, 73, 61, 69, 65, 76, 59, 69, 63, 55, 65, 68, 69, 70, 62, 68, 76, 63, 72, 64, 
65, 64, 74, 77, 67, 60, 74, 72, 72, 74, 64, 59, 70, 63, 61, 72, 76, 61, 70, 68, 68, 75, 66, 64, 78, 72, 63, 72, 72, 69, 70, 71, 74, 62, 69, 73, 72, 69, 63, 80, 70, 63, 62, 73, 56, 64, 64, 67, 72, 74, 58, 65, 77, 70, 73, 70, 
68, 70, 78, 86, 76, 77, 64, 65, 67, 78, 67, 67, 71, 61, 62, 66, 76, 64, 75, 78, 69, 66, 65, 69, 71, 60, 65, 72, 79, 76, 70, 58, 66, 74, 70, 68, 63, 72, 63, 66, 69, 78, 78, 76, 60, 77, 75, 68, 67, 73, 70, 72, 69, 70, 68, 64]

mean_exp = 6

# Please run %matplotlib qt before %run StatConnProject.py

ca3 = CA3_Model(1008)
routine(ca3,[0.75,0.25],0.7,shuff_counts,6)
	
        
