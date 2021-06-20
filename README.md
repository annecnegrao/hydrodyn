# hydrodyn
1D Hydrodynamic Model for urban channels that I developed during my master's degree.

Developed by _Anne Caroline Negrão_
**annecnegrao@gmail.com**

--------

## _AkanMelinaLPI-Gregorio.f90_

Simulation of Gregorio creek considering the following channel bottoms:
- simple bottom (fundo_simples)
- bottom with steps (fundo_degrau)

Into this directory you will find the main information about the channel:
fundo_simnples.txt
fundo_degraus.txt 
Where x is the horizontal distance in meters, S0 is the slope and
z is the depth.

Into info_canal.txt there is some information for the symulation.

Into arc_eventos.txt there is a list with the file name of events that
will be simulated.

If heating process (1) is activate then:
    Generation of initial conditions from a linear waterline
        create directories for each event
        starts a linear waterline from the level of the first instant of the event
        heating the initial condition

If simulation process (2) is activate then:
    Simulation of Gregorio creek with heated initial condition
    Calculation of the Nash-Sutcliffe coefficient

## Inputs:
    info_canal.txt
    fundo_simplificado.txt or fundo_real.txt
    arq_eventos.txt or arq_eventos_cut.txt
    CI_aquecido.txts

## Outputs:
    arq_eventos_cut.txt
    ./<event_directory>/h_y_obs.txt
    ./<event_directory>/h_y_obs.txt
    ./<event_directory>/CI_linear.txt
    ./<event_directory>/CI_aquecido.txt
    ./<event_directory>/linha.gnu
    ENS.txt
    ./<event_directory>/res_xh.txt
    ./<event_directory>/res_xQ.txt
    ./<event_directory>/res_xy.txt
    ./<event_directory>/res_tQ.txt
    ./<event_directory>/res_ty.txt


    