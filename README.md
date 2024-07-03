/** Processor calculates MC decay angles from W boson.
 *  PROTOCOL
 * 1) Goes through MC Particles
 * 
 * 2) [ON SHELL] W has PDG abs(24) ---> Get four-momentum && Charge
 *  -  Production angle of W- (CosineTheta)
 *  -  LEPTONIC and HADRONIC decay discrimination (Directly from W daughters).
 *  -  LEPTONIC DECAY FROM ---> W- or W+ (cross check by counting) ---> For decay angles
 *  -  HADRONIC DECAY FROM ---> W- or W+ (cross check by counting)
 *  -  DECAY ANGLES ---> Use starVector ---> *theta and *phi
 * 
 * 3) [OFF SHELL] W reconstructed from it's daughters ---> Get four-momentum && Charge
 *  -  Production angle of W- (CosineTheta)
 *  -  LEPTONIC and HADRONIC decay discrimination (pdg > or < than 10) 
 *  -  LEPTONIC FROM ---> W- or W+ (cross check by counting) ---> Use for decay angles
 *  -  HADRONIC FROM ---> W- or W+ (cross check by counting)
 *  -  DECAY ANGLES ---> Use starVector -> *theta and *phi
 * 
 *  @author A.Silva, DESY
 */
   @author A.Silva, DESY
