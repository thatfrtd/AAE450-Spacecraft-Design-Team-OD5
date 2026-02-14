function [delta_v, t] = Rendezvouz_and_Docking()  %total delta V used for rendezvus and docking

   arguments (Input)

   end
    
    delta_v = 20 * 10 ^ -3;  %[km/s] estimate from Richard et. al, 2013

    t = 6 * 3600;  %[s] Rendezvous Time

end