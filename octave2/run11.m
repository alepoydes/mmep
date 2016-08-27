more off
system=latticetriagonal([31,31,1],0.1,0.1,1,0.3);
fd=project(randn([system.size,3]));
plotfield(system,fd,'-k');
topcharge(fd)/4/pi
topcharge2(fd)/4/pi