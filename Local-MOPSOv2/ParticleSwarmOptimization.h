#ifndef CONTINUOUSPARTICLESWARMOPT_H
#define CONTINUOUSPARTICLESWARMOPT_H

#include "Swarm.h"

#define NOTCHANGED      -1

/* -10<=V<=10 */
#define sigmoid(V,M)            M/(1.0 + exp(-V))

/*Population, EvaluateParticle, TerminationCriteria*/
void SwarmOptimization(Swarm *, void (*)(Particle *), int (*)(Swarm *));

//Particle *ContinuousBestGlobalParticle(Population *);
#endif
