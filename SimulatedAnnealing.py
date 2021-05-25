#  SA Implementation for TSP
#  Author: Konrad Eilers

from random import randint
from random import random
from random import shuffle
from random import sample
from random import uniform
import math
import os
import numpy as np

import matplotlib.pyplot as plt


class SA:

    def __init__(self, x0, step_max, t_min, t_max, cooling_schedule_type, alpha, tau):
        
        # Starting conditions
        self.t = t_max # current time 
        self.t_max = t_max
        self.t_min = t_min
        
        self.step_max = step_max
        self.step_max = step_max

        self.hist = [] # store of history

        self.current_state = x0
        self.current_energy = self.cost_function(x0)
        self.average_energy = 0
        self.best_state = self.current_state
        self.best_energy = self.current_energy
        self.alpha = alpha
        self.tau = tau

        if cooling_schedule_type not in ['lin_mult', 'lin', 'geometric', 'quad_mult']:
            raise NotImplementedError
        else:
            self.type = cooling_schedule_type

    def safe_exp(self, x):
        try:
            return math.exp(x)
        except:
            return 0
    
    def cooling_schedule(self, step):

        if self.type == 'lin_mult':
            return self.t_max /  (1 + self.alpha * step)

        if self.type == 'lin':
            return self.t_max - self.alpha * step
 
        if self.type == 'geometric':
            return self.t_max * self.alpha ** step

        if self.type == 'quad_mult':
            return self.t_max /  (1 + self.alpha * step**2)

    def get_neighbor(self):
        """ 
        Return neighboring solution to current point by swapping two points
        """

        p0 = randint(0, len(self.current_state)-1)
        p1 = randint(0, len(self.current_state)-1)

        neighbor = self.current_state[:]
        neighbor[p0], neighbor[p1] = neighbor[p1], neighbor[p0]

        return neighbor

    def calc_euclidean(self, p1, p2):    
        return ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**0.5
    
    def cost_function(self, nodes):
            dist = 0
            for i in range(len(nodes)-1): 
                dist += self.calc_euclidean(nodes[i+1], nodes[i])
            dist += self.calc_euclidean(nodes[0], nodes[-1])
            return dist

    def optimization(self):
        """
        # Implementaiton of optimization routine
        """

        self.step, self.accept = 1, 0

        while self.step < self.step_max and self.t >= self.t_min:

            energy_store = []

            for tau in range(self.tau):

                # get neighbor
                proposed_neighbor = self.get_neighbor() 

                # check energy level of neighbor
                E_n = self.cost_function(proposed_neighbor)
                dE = E_n - self.current_energy

                # determine if we should accept the current neighbor
                if random() < self.safe_exp(-max(0,dE) / self.t): 
                    self.current_energy = E_n
                    self.current_state = proposed_neighbor[:]
                    self.accept += 1

                    energy_store.append(self.current_energy)
                
                # check if the current neighbor is best solution so far
                if E_n < self.best_energy:
                    self.best_energy = E_n
                    self.best_state = proposed_neighbor[:]
            
            self.average_energy = np.mean(energy_store)

            # record keeping
            self.hist.append([
                self.step,
                self.t,
                self.average_energy,
                self.best_energy])
            
            # below should be outside the tau loop 

            self.t = self.cooling_schedule(self.step)
            self.step += 1
        
        return self.hist
    

class TravellingSalesman(SA):

    def __init__(self, x0, step_max, t_min, t_max, cooling_schedule_type, alpha, tau):
        super().__init__(x0, step_max, t_min, t_max, cooling_schedule_type, alpha, tau)
        
    
    def run(self):
        self.optimization()
        return self.hist, self.best_state

def data_generator_cirlce():
    n_pts = 50
    dr = (2 * math.pi) / n_pts

    x0 = []
    for i in range(n_pts):
        radians = dr * i  
        x0.append([math.cos(radians), math.sin(radians)])
    
    return sample(x0, n_pts)

def random_data_generator():
    n_pts = 50

    x0 = []
    for i in range(n_pts):
        x0.append([uniform(-1,1), uniform(-1,1)])
    
    return sample(x0, n_pts)


def main():

    # build directory for plots
    directory = 'simulated_anealing_plots'
    if not os.path.exists(directory):
        os.makedirs(directory)

    data_points = data_generator_cirlce()
    history, best_state = TravellingSalesman(x0=data_points,
                                            step_max=1000, 
                                            t_min=0.01, 
                                            t_max=20, 
                                            cooling_schedule_type='quad_mult',  
                                            alpha=0.001, 
                                            tau = 1000).run()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize = (20,4))

    data_points = data_points + [data_points[0]]
    best_state = best_state + [best_state[0]]

    fig.suptitle('Traveling Salesman Problem')
    ax1.plot(*zip(*data_points))
    ax1.set_title('Initial Situation')
    ax2.plot(*zip(*best_state))
    ax2.set_title('After Opimization')

    plt.savefig('simulated_anealing_plots/before_after_opt.png')

    temp = []
    step = []
    average_state = []
    best_state = []
    for i in range(len(history)):
        step.append(history[i][0])
        temp.append(history[i][1])
        average_state.append(history[i][2])
        best_state.append(history[i][3])

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, figsize = (10,8))
    fig.suptitle('Temperature, Average State, Best State')
    ax1.plot(step, temp)
    ax1.set_ylabel('Temperature')
    ax2.plot(step, average_state)
    ax2.set_ylabel('Energy')
    ax3.plot(step, best_state)
    ax3.set_ylabel('Energy')
    ax3.set_xlabel('Steps')

    plt.savefig('simulated_anealing_plots/step_history.png')

if __name__ == "__main__":
    main()

