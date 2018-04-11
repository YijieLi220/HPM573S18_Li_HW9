import scr.RandomVariantGenerators as rndClasses
import scr.StatisticalClasses as Stat
import scr.SamplePathClasses as PathCls
from enum import Enum
import scr.FigureSupport as Figs

TRANS_MATRIX = [
    [0.75,  0.15,   0,     0.1],   # WELL
    [0,     0,      1,     0],   # STROKE
    [0,     0.25,   0.55,  0.2],   # P_STROKE
    ]

TRANS_MATRIX_THERAPY = [
    [0.75,  0.15,    0,      0.1],   # WELL
    [0,     0,       1,      0],   # STROKE
    [0,     0.1625,  0.701,  0.1365],   # P_STROKE
    ]

TRANS=[[
    [0.75,  0.15,    0,     0.1],   # WELL
    [0,     0,       1,     0],   # STROKE
    [0,     0.25,    0.55,  0.2],   # P_STROKE
    ], [
    [0.75,  0.15,    0,      0.1],   # WELL
    [0,     0,       1,      0],   # STROKE
    [0,     0.1625,  0.701,  0.1365],   # P_STROKE
    ]
]


class HealthStats:
    """ health states of patients with risk of stroke """
    WELL = 0
    STROKE = 1
    P_STROKE = 2
    DEATH = 3

class THERAPY_OR_NOT (Enum):
    WITHOUT=0
    WITH=1


class Patient:
    def __init__(self, id, THERAPY):
        """ initiates a patient
        :param id: ID of the patient
        :param parameters: parameter object
        """

        self._id = id
        # random number generator for this patient
        self._rng = None
        self.healthstat=0
        self.survival=0
        self.THERA = THERAPY
        self.STROKE=0

    def simulate(self, sim_length):
        """ simulate the patient over the specified simulation length """

        # random number generator for this patient
        self._rng = rndClasses.RNG(self._id)

        k = 0  # current time step

        # while the patient is alive and simulation length is not yet reached
        while self.healthstat!=3 and k  < sim_length:
            # find the transition probabilities of the future states
            trans_probs = TRANS[self.THERA][self.healthstat]
            # create an empirical distribution
            empirical_dist = rndClasses.Empirical(trans_probs)
            # sample from the empirical distribution to get a new state
            # (returns an integer from {0, 1, 2, ...})
            new_state_index = empirical_dist.sample(self._rng)
            if self.healthstat==1:
                self.STROKE+=1
            # update health state
            self.healthstat =new_state_index[0]
            # increment time step
            k += 1
        self.survival=k+1

    def get_survival_time(self):
        """ returns the patient's survival time"""
        return self.survival

    def get_time_to_STROKE(self):
        """ returns the patient's survival time"""
        return self.STROKE

##for question#3
class Cohort():
    def __init__(self,id,THERAPY):
        self._initial_pop_size=2000
        self.survivaltime=[]
        self.id=id
        self.THERAPY=THERAPY
        self.STROKE=[]

    def simulate(self):
        for i in range(self._initial_pop_size):
            patient=Patient(self.id*self._initial_pop_size+i,self.THERAPY)
            patient.simulate(1000)
            self.survivaltime.append(patient.get_survival_time())
            self.STROKE.append(patient.get_time_to_STROKE())

    def get_survival_time(self):
        return self.survivaltime

    def get_time_to_STROKE(self):
        """ returns the patient's survival time"""
        return self.STROKE


class CohortOutputs:
    def __init__(self, simulated_cohort):
        """ extracts outcomes of a simulated cohort
        :param simulated_cohort: a cohort after being simulated"""

        self._simulatedCohort = simulated_cohort

    def get_ave_survival_time(self):
        """ returns the average survival time of patients in this cohort """
        return sum(self._simulatedCohort.get_survival_time()) / len(self._simulatedCohort.get_survival_time())

    def get_survival_curve(self):
        """ returns the sample path for the number of living patients over time """

        # find the initial population size
        n_pop = 2000
        # sample path (number of alive patients over time)
        n_living_patients = PathCls.SamplePathBatchUpdate('# of living patients', 0, n_pop)

        # record the times of deaths
        for obs in self._simulatedCohort.get_survival_time():
            n_living_patients.record(time=obs, increment=-1)

        return n_living_patients

    def get_survival_times(self):
        """ :returns the survival times of the patients in this cohort"""
        return self._simulatedCohort.get_survival_time()



cohort_1=Cohort(1,THERAPY_OR_NOT.WITHOUT.value)
cohort_1.simulate()

cohort_2=Cohort(1,THERAPY_OR_NOT.WITH.value)
cohort_2.simulate()

sum_stat = Stat.SummaryStat("dsa",cohort_1.get_survival_time())
ExpectedCI=sum_stat.get_t_CI(0.05)
meansurvival=sum_stat.get_mean()
print(meansurvival,ExpectedCI)

sum_stat_1 = Stat.SummaryStat("dsa",cohort_2.get_survival_time())
ExpectedCI_1=sum_stat_1.get_t_CI(0.05)
meansurvival_1=sum_stat_1.get_mean()
print(meansurvival_1,ExpectedCI_1)
##PLOT:
cohortOUT_1=CohortOutputs(cohort_1)
cohortOUT_2=CohortOutputs(cohort_2)
survival_curves = [
    cohortOUT_1.get_survival_curve(),
    cohortOUT_2.get_survival_curve()
]

survival_curves = [
    cohortOUT_1.get_survival_curve(),
    cohortOUT_2.get_survival_curve()
]

PathCls.graph_sample_paths(
    sample_paths=survival_curves,
    title='Survival curve',
    x_label='Simulation time step',
    y_label='Number of alive patients',
    legends=['No Drug', 'With Drug']
)

set_STROKE_times = [
    cohort_1.get_time_to_STROKE(),
    cohort_2.get_time_to_STROKE()
]

Figs.graph_histograms(
    data_sets=set_STROKE_times,
    title='Histogram of patient STROKE time',
    x_label='Survival time',
    y_label='Counts',
    bin_width=1,
    legend=['WELL', 'POST DRUG'],
    transparency=0.6
)