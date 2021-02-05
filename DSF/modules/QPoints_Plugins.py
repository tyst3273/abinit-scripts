"""
This file should be plugins to generate lists of Qpoints. Basically all that is required is to create the 
params.Qpoints attribute that later is passed to Phonopy or anaddb

params.Qpoints should be a (num_Qpoints x 3) numpy array where num_Qpoints is the number of Qpoints in the BZ path

User should call compute_reduced_qpoints at the end of the method 

"""

import numpy as np

class plugins:



    def brillouin_zone_slice(self,params,Q_path,num_Qpoints):

        params.Qpoints = np.zeros((num_Qpoints,3))
        params.Qpoints[:,0] = np.linspace(Q_path[0][0],Q_path[1][0],num_Qpoints)
        params.Qpoints[:,1] = np.linspace(Q_path[0][1],Q_path[1][1],num_Qpoints)
        params.Qpoints[:,2] = np.linspace(Q_path[0][2],Q_path[1][2],num_Qpoints)
        params.num_Qpoints = params.Qpoints.shape[0]

        self.compute_reduced_qpoints(params)



    def brillouin_zone_path(self,params,Q_list=[[[0.0,0.0,0.0],[0.0,0.5,0.5]],[[0.5,0.5,0.5],[0.0,0.0,0.0]]],
                                            num_Qpoints_in_shortest_segment=101):
        """
        give endpoints of a slice of the 1st BZ. can have disjoint segments together by giving lists of lists 
        of endpoints
        """

        num_segments = len(Q_list)
        split_segments = False
        for l in Q_list:
            if type(l[0]) == list:
                split_segments = True
                break

        Qpoints = []
        Qlengths = []
        min_len = 1e100

        # kinda a stupid amount of repeated loops, but its not really expensive
        if split_segments == True:
            for l in Q_list:
                lengths = []
                num_endpoints = len(l)
                for i in range(num_endpoints-1):
                    lengths.append(np.sqrt((l[i][0]-l[i+1][0])**2+(l[i][1]-l[i+1][1])**2+(l[i][2]-l[i+1][2])**2))
                if min(lengths[:]) <= min_len:
                    min_len = min(lengths[:])
                Qlengths.append(lengths)

            for i in range(len(Qlengths)):
                Qlengths[i][:] = Qlengths[i][:]/min_len*num_Qpoints_in_shortest_segment
                Qlengths[i] = [int(j) for j in Qlengths[i]]

            for i in range(len(Qlengths)):
                points = []
                for j in range(len(Qlengths[i])):
                    tmp = np.zeros((Qlengths[i][j],3))
                    if j == len(Qlengths[i])-1:
                        tmp[:,0] = np.linspace(Q_list[i][j][0],Q_list[i][j+1][0],Qlengths[i][j])
                        tmp[:,1] = np.linspace(Q_list[i][j][1],Q_list[i][j+1][1],Qlengths[i][j])
                        tmp[:,2] = np.linspace(Q_list[i][j][2],Q_list[i][j+1][2],Qlengths[i][j])
                    else:
                        tmp[:,0] = np.linspace(Q_list[i][j][0],Q_list[i][j+1][0],Qlengths[i][j]+1)[:Qlengths[i][j]]
                        tmp[:,1] = np.linspace(Q_list[i][j][1],Q_list[i][j+1][1],Qlengths[i][j]+1)[:Qlengths[i][j]]
                        tmp[:,2] = np.linspace(Q_list[i][j][2],Q_list[i][j+1][2],Qlengths[i][j]+1)[:Qlengths[i][j]]
                    points.append(tmp)
                Qpoints.append(points)

            total_num_Qpoints = 0
            for l in Qlengths:
                for j in l:
                    total_num_Qpoints = total_num_Qpoints+j
            params.Qpoints = np.zeros((total_num_Qpoints,3))
            count = 0
            for i in range(len(Qlengths)):
                for j in range(len(Qlengths[i])):
                    params.Qpoints[count:count+Qlengths[i][j],:] = Qpoints[i][j]
                    count = count+Qlengths[i][j]

        else:
            lengths = []
            for i in range(len(Q_list)-1):
                lengths.append(np.sqrt((Q_list[i][0]-Q_list[i+1][0])**2+
                    (Q_list[i][1]-Q_list[i+1][1])**2+(Q_list[i][2]-Q_list[i+1][2])**2))
            if min(lengths[:]) <= min_len:
                    min_len = min(lengths[:])
            Qlengths.extend(lengths)
            Qlengths[:] = Qlengths[:]/min_len*num_Qpoints_in_shortest_segment
            Qlengths = [int(j) for j in Qlengths]

            total_num_Qpoints = 0
            for l in Qlengths:
                total_num_Qpoints = total_num_Qpoints+l

            params.Qpoints = np.zeros((total_num_Qpoints,3))
            count = 0
            for i in range(len(Qlengths)):
                if i == len(Qlengths)-1:
                    params.Qpoints[count:count+Qlengths[i],0] = np.linspace(Q_list[i][0],Q_list[i+1][0],Qlengths[i])
                    params.Qpoints[count:count+Qlengths[i],1] = np.linspace(Q_list[i][1],Q_list[i+1][1],Qlengths[i])
                    params.Qpoints[count:count+Qlengths[i],2] = np.linspace(Q_list[i][2],Q_list[i+1][2],Qlengths[i])
                else:
                    params.Qpoints[count:count+Qlengths[i],0] = np.linspace(Q_list[i][0],Q_list[i+1][0],
                            Qlengths[i]+1)[:Qlengths[i]]
                    params.Qpoints[count:count+Qlengths[i],1] = np.linspace(Q_list[i][1],Q_list[i+1][1],
                            Qlengths[i]+1)[:Qlengths[i]]
                    params.Qpoints[count:count+Qlengths[i],2] = np.linspace(Q_list[i][2],Q_list[i+1][2],
                            Qlengths[i]+1)[:Qlengths[i]]
                count = count+Qlengths[i]

#        np.savetxt('QPoints_list',params.Qpoints,fmt='%3.9f')        
        params.num_Qpoints = params.Qpoints.shape[0]
        self.compute_reduced_qpoints(params)



    def explicit_list(self,params,Qpoints_list):

        params.num_Qpoints = len(Qpoints_list)
        params.Qpoints = np.zeros((params.num_Qpoints,3))
        for i in range(params.num_Qpoints):
            params.Qpoints[i,:] = Qpoints_list[i]

        self.compute_reduced_qpoints(params)



    def read_file(self,params,filename):
    
        params.Qpoints = np.loadtxt(filename)
        params.num_Qpoints = params.Qpoints.shape[0]

        self.compute_reduced_qpoints(params)



######################################################################################

    def compute_reduced_qpoints(self,params):

        params.tau = np.round(params.Qpoints)
        params.qpoints_plus = params.Qpoints-params.tau
        params.qpoints_minus = params.tau-params.Qpoints
        





















