import datetime
import logging

import numpy as np
import scipy
import tqdm
from petsc4py import PETSc

from ._Muskingum import Muskingum


class MuskingumNumerical(Muskingum):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _route(self,
               dates: np.array,
               runoffs: np.array,
               q_t: np.array,
               r_t: np.array,
               inflow_t: np.array) -> np.array:
        self.make_lhs_matrix()
        outflow_array = np.zeros((self.num_outflow_steps, self.A.shape[0]))

        logging.info('Creating PETSc arrays')
        self.lhs = scipy.sparse.csr_matrix(self.lhs)
        A = PETSc.Mat().createAIJ(size=self.lhs.shape, csr=(self.lhs.indptr, self.lhs.indices, self.lhs.data))
        x = PETSc.Vec().createSeq(size=self.lhs.shape[0])
        b = PETSc.Vec().createSeq(size=self.lhs.shape[0])

        # Define a KSP (Krylov Subspace Projection) solver
        ksp = PETSc.KSP().create()
        ksp.setType(self.conf.get('petsc_ksp_type', 'richardson'))
        ksp.setTolerances(rtol=1e-5)
        ksp.setOperators(A)

        # Define a preconditioner if specified in config
        if self.conf.get('petsc_pc_type', ''):
            pc = PETSc.PC().create()
            pc.setType(self.conf['petsc_pc_type'])
            pc.setOperators(A)
            pc.setUp()
            ksp.setPC(pc)

        logging.info('Performing routing solver iterations')
        t1 = datetime.datetime.now()
        for inflow_time_step, inflow_end_date in enumerate(tqdm.tqdm(dates, desc='Runoff Routed')):
            r_t = runoffs[inflow_time_step, :]
            interval_flows = np.zeros((self.num_routing_substeps_per_outflow, self.A.shape[0]))
            for routing_substep_iteration in range(self.num_routing_substeps_per_outflow):
                inflow_tnext = (self.A @ q_t) + r_t
                rhs = (self.c1 * inflow_t) + (self.c2 * r_t) + (self.c3 * q_t)
                inflow_t = inflow_tnext
                b.setArray(rhs)
                ksp.solve(b, x)
                q_t = x.getArray()
                interval_flows[routing_substep_iteration, :] = q_t
            interval_flows = np.mean(interval_flows, axis=0)
            interval_flows = np.round(interval_flows, decimals=2)
            outflow_array[inflow_time_step, :] = interval_flows
        t2 = datetime.datetime.now()
        logging.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')

        self.qinit = q_t
        self.rinit = r_t

        logging.info('Cleaning up PETSc objects')
        A.destroy()
        x.destroy()
        b.destroy()
        ksp.destroy()
        pc.destroy() if self.conf.get('petsc_pc_type', '') else None

        return outflow_array
