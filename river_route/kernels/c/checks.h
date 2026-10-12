/*
 * What routing.c calls in checks.c before routing a call. It is hidden from the native library's exports, which are
 * only those of routing.h.
 */
#ifndef RIVER_ROUTE_CHECKS_H
#define RIVER_ROUTE_CHECKS_H

#include "routing.h"

/* Check every river of a call and start it: each river's discharge row at the discharge of its last sub-reach. */
__attribute__((visibility("hidden"))) int start_job(const struct rr_network *network_pointer,
                                                     const struct rr_runoff *runoff_pointer,
                                                     const struct rr_job *job_pointer);

#endif
