import datetime
import logging

import numpy as np
import tqdm

from ._Muskingum import Muskingum


class MuskingumAnalytical(Muskingum):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _route(self,
               dates: np.array,
               runoffs: np.array,
               q_t: np.array,
               r_t: np.array,
               inflow_t: np.array) -> np.array:
        outflow_array = np.zeros((self.num_outflow_steps, self.A.shape[0]))

        logging.info('Performing routing computation iterations')
        t1 = datetime.datetime.now()
        for inflow_time_step, inflow_end_date in enumerate(tqdm.tqdm(dates, desc='Runoff Routed')):
            r_t = runoffs[inflow_time_step, :]
            interval_flows = np.zeros((self.num_routing_substeps_per_outflow, self.A.shape[0]))
            for routing_substep_iteration in range(self.num_routing_substeps_per_outflow):
                inflow_tnext = (self.A @ q_t) + r_t
                q_t = self.lhsinv @ ((self.c1 * inflow_t) + (self.c2 * r_t) + (self.c3 * q_t))
                interval_flows[routing_substep_iteration, :] = q_t
                inflow_t = inflow_tnext
            interval_flows = np.mean(interval_flows, axis=0)
            interval_flows = np.round(interval_flows, decimals=2)
            outflow_array[inflow_time_step, :] = interval_flows
        t2 = datetime.datetime.now()
        logging.info(f'Routing completed in {(t2 - t1).total_seconds()} seconds')

        self.qinit = q_t
        self.rinit = r_t
        return outflow_array
