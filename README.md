/** Processor calculates decay angles from W boson.
 *  PROTOCOL
 * 1) Goes through MC Particles
 * 
 * 2) ON SHELL W has PDG abs(24) -> production angle of leptonic W 
 *  -  Discrimination for LEPTONIC and HADRONIC decay -> GetDaughters() -> GetPDG 
 *  -  If this a lepton PDG, this is a leptonic W -> production angle of leptonic Ws
 *  -  Use starVector -> *theta and *phi
 * 
 * 3) OFF SHELL W reconstructed from it's daughters:
 *  3.1) Discrimination for LEPTONIC and HADRONIC decay (pdg > or < than 10) -> production angle of leptonic Ws
 *  3.2) LEPTONIC FROM ---> W- or W+ (cross check by counting)
 *  3.3) HADRONIC FROM ---> W- or W+ (cross check by counting)
 *  3.4) Use starVector -> *theta and *phi
 * 
 *  @author A.Silva, DESY
 */
   @author A.Silva, DESY
