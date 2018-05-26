# Makefile created by mkmf 19.2.0


.DEFAULT:
	-echo $@ does not exist.
all: a.out
F77KINDS.o: ./F77KINDS.F90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FPPFLAGS) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./F77KINDS.F90
acoustics.o: ./acoustics.f90 fishy.o geomc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./acoustics.f90
aerate.o: ./aerate.f90 global.o main.o kinetic.o trans.o screenc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./aerate.f90
balances.o: ./balances.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./balances.f90
barea.o: ./barea.f90 F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./barea.f90
bioenergetics.o: ./bioenergetics.f90 F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./bioenergetics.f90
calcfric.o: ./calcfric.f90 global.o geomc.o eddy.o kinetic.o logicc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./calcfric.f90
calculate_az.o: ./calculate_az.f90 geomc.o global.o trans.o eddy.o kinetic.o macrophytec.o logicc.o main.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./calculate_az.f90
cdfunc.o: ./cdfunc.f90 structures.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cdfunc.f90
#ce_qual_w2.o: ./ce_qual_w2.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o cemavars.o initialvelocity.o envirpmod.o bioenergetics.o F77KINDS.o
#	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ce_qual_w2.f90
cema_w2_input.o: ./cema_w2_input.f90 main.o global.o geomc.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cema_w2_input.f90
cemabubblesrelease.o: ./cemabubblesrelease.f90 main.o global.o geomc.o cemavars.o screenc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemabubblesrelease.f90
cemabubblesreleaseturbulence.o: ./cemabubblesreleaseturbulence.f90 main.o global.o geomc.o screenc.o kinetic.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemabubblesreleaseturbulence.f90
cemabubblesrisevelocity.o: ./cemabubblesrisevelocity.f90 F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemabubblesrisevelocity.f90
cemabubblestransport.o: ./cemabubblestransport.f90 main.o global.o geomc.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemabubblestransport.f90
cemabubblesturbulence.o: ./cemabubblesturbulence.f90 main.o global.o geomc.o screenc.o kinetic.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemabubblesturbulence.f90
cemabubbwattransfer.o: ./cemabubbwattransfer.f90 main.o global.o geomc.o screenc.o kinetic.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemabubbwattransfer.f90
cemacalculaterisevelocity.o: ./cemacalculaterisevelocity.f90 main.o global.o geomc.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemacalculaterisevelocity.f90
cemacomputeturbidity.o: ./cemacomputeturbidity.f90 global.o screenc.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemacomputeturbidity.f90
cemadisgasphasedistribution.o: ./cemadisgasphasedistribution.f90 cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemadisgasphasedistribution.f90
cemafftlayercode.o: ./cemafftlayercode.f90 main.o global.o geomc.o screenc.o rstart.o eddy.o logicc.o tvdc.o kinetic.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemafftlayercode.f90
cemaoutputroutines.o: ./cemaoutputroutines.f90 main.o global.o screenc.o geomc.o kinetic.o namesc.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemaoutputroutines.f90
cemasedimentdiagenesis.o: ./cemasedimentdiagenesis.f90 main.o global.o geomc.o screenc.o rstart.o eddy.o logicc.o tvdc.o kinetic.o surfhe.o cemavars.o trans.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemasedimentdiagenesis.f90
cemasedimentmodelw2.o: ./cemasedimentmodelw2.f90 main.o global.o geomc.o screenc.o rstart.o eddy.o logicc.o tvdc.o kinetic.o trans.o cemavars.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemasedimentmodelw2.f90
cemavars.o: ./cemavars.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./cemavars.f90
deallocate_cema.o: ./deallocate_cema.f90 cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./deallocate_cema.f90
deg_c.o: ./deg_c.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./deg_c.f90
deg_f.o: ./deg_f.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./deg_f.f90
density.o: ./density.f90 logicc.o global.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./density.f90
depthcrit.o: ./depthcrit.f90 structures.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./depthcrit.f90
eddy.o: ./eddy.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./eddy.f90
enable.o: ./enable.f90 
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./enable.f90
endsimulation.o: ./endsimulation.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o initialvelocity.o bioenergetics.o tridiag_v.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./endsimulation.f90
envirp.o: ./envirp.f90 global.o main.o namesc.o screenc.o tvdc.o rstart.o geomc.o envirpmod.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./envirp.f90
envirpmod.o: ./envirpmod.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./envirpmod.f90
fimpbr.o: ./fimpbr.f90 fishy.o geomc.o global.o surfhe.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./fimpbr.f90
findbranch.o: ./findbranch.f90 fishy.o global.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./findbranch.f90
findnewbr.o: ./findnewbr.f90 fishy.o geomc.o global.o screenc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./findnewbr.f90
fish.o: ./fish.f90 surfhe.o fishy.o gdayc.o screenc.o geomc.o global.o main.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./fish.f90
fishhabitat.o: ./fishhabitat.f90 global.o main.o screenc.o kinetic.o tvdc.o namesc.o logicc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./fishhabitat.f90
fishoutput.o: ./fishoutput.f90 fishy.o global.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./fishoutput.f90
fishplot.o: ./fishplot.f90 fishy.o geomc.o global.o gdayc.o screenc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./fishplot.f90
fishy.o: ./fishy.f90 F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./fishy.f90
gas_transfer.o: ./gas_transfer.f90 global.o geomc.o kinetic.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./gas_transfer.f90
gasbubblesformation.o: ./gasbubblesformation.f90 global.o screenc.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./gasbubblesformation.f90
gasdev.o: ./gasdev.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./gasdev.f90
gate_flow.o: ./gate_flow.f90 structures.o global.o geomc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./gate_flow.f90
gdayc.o: ./gdayc.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./gdayc.f90
geomc.o: ./geomc.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./geomc.f90
global.o: ./global.f90 F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./global.f90
gregorian_date.o: ./gregorian_date.f90 gdayc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./gregorian_date.f90
grid_area1.o: ./grid_area1.f90 global.o geomc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./grid_area1.f90
grid_area2.o: ./grid_area2.f90 global.o geomc.o rstart.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./grid_area2.f90
gridplot.o: ./gridplot.f90 fishy.o global.o geomc.o screenc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./gridplot.f90
heat_exchange.o: ./heat_exchange.f90 global.o gdayc.o surfhe.o tvdc.o shadec.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./heat_exchange.f90
hydroinout.o: ./hydroinout.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./hydroinout.f90
init.o: ./init.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o initialvelocity.o bioenergetics.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./init.f90
initcond.o: ./initcond.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./initcond.f90
initgeom.o: ./initgeom.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o porosityc.o macrophytec.o rstart.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./initgeom.f90
initial_u_velocity.o: ./initial_u_velocity.f90 global.o geomc.o initialvelocity.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./initial_u_velocity.f90
initial_water_level.o: ./initial_water_level.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o initialvelocity.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./initial_water_level.f90
initializebedconsolidationfiles.o: ./initializebedconsolidationfiles.f90 main.o global.o screenc.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./initializebedconsolidationfiles.f90
initialvelocity.o: ./initialvelocity.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./initialvelocity.f90
input.o: ./input.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o initialvelocity.o bioenergetics.o tridiag_v.o 
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./input.f90
interconst.o: ./interconst.f90 fishy.o global.o geomc.o screenc.o kinetic.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./interconst.f90
interflowf.o: ./interflowf.f90 fishy.o global.o geomc.o screenc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./interflowf.f90
kinetic.o: ./kinetic.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./kinetic.f90
kinetics.o: ./kinetics.f90 screenc.o global.o kinetic.o geomc.o tvdc.o logicc.o surfhe.o macrophytec.o zooplanktonc.o main.o cemavars.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./kinetics.f90
lateral_velocity.o: ./lateral_velocity.f90 fishy.o global.o geomc.o main.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./lateral_velocity.f90
layeraddsub.o: ./layeraddsub.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./layeraddsub.f90
lin_sys.o: ./lin_sys.f90 F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./lin_sys.f90
logicc.o: ./logicc.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./logicc.f90
lubksb.o: ./lubksb.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./lubksb.f90
ludcmp.o: ./ludcmp.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./ludcmp.f90
macrophyte_friction.o: ./macrophyte_friction.f90 geomc.o global.o macrophytec.o porosityc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./macrophyte_friction.f90
macrophytec.o: ./macrophytec.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./macrophytec.f90
main.o: ./main.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./main.f90
mannings_eqn.o: ./mannings_eqn.f90 global.o geomc.o eddy.o logicc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./mannings_eqn.f90
#msclib.o: ./msclib.f90
#	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./msclib.f90
namesc.o: ./namesc.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./namesc.f90
normal_depth.o: ./normal_depth.f90 global.o geomc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./normal_depth.f90
open_channel_initialize.o: ./open_channel_initialize.f90 global.o structures.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./open_channel_initialize.f90
output.o: ./output.f90 global.o gdayc.o geomc.o kinetic.o tvdc.o namesc.o logicc.o macrophytec.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./output.f90
outputa.o: ./outputa.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o bioenergetics.o cemavars.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./outputa.f90
part_transport.o: ./part_transport.f90 fishy.o trans.o global.o geomc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./part_transport.f90
pipe_flow_initialize.o: ./pipe_flow_initialize.f90 global.o geomc.o structures.o screenc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./pipe_flow_initialize.f90
porosity.o: ./porosity.f90 geomc.o global.o macrophytec.o porosityc.o logicc.o screenc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./porosity.f90
porosityc.o: ./porosityc.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./porosityc.f90
random.o: ./random.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./random.f90
read_fish_data.o: ./read_fish_data.f90 fishy.o screenc.o main.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./read_fish_data.f90
readbedconsolidationfiles.o: ./readbedconsolidationfiles.f90 main.o global.o screenc.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./readbedconsolidationfiles.f90
restart_output.o: ./restart_output.f90 global.o screenc.o rstart.o gdayc.o geomc.o kinetic.o tvdc.o zooplanktonc.o eddy.o main.o envirpmod.o logicc.o structures.o macrophytec.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./restart_output.f90
rstart.o: ./rstart.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./rstart.f90
rtbis.o: ./rtbis.f90 global.o eddy.o screenc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./rtbis.f90
run_w2.o: ./run_w2.f90 global.o main.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./run_w2.f90
screen_update.o: ./screen_update.f90 geomc.o global.o gdayc.o screenc.o surfhe.o tvdc.o logicc.o namesc.o structures.o main.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./screen_update.f90
screenc.o: ./screenc.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./screenc.f90
selective.o: ./selective.f90 selective1.o main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./selective.f90
selective1.o: ./selective1.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./selective1.f90
selective1usgs.o: ./selective1usgs.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./selective1usgs.f90
selectiveinit.o: ./selectiveinit.f90 selective1.o main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./selectiveinit.f90
selectiveinitusgs.o: ./selectiveinitusgs.f90 selective1usgs.o main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./selectiveinitusgs.f90
selectiveusgs.o: ./selectiveusgs.f90 selective1usgs.o main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./selectiveusgs.f90
selwc.o: ./selwc.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./selwc.f90
semilog.o: ./semilog.f90 eddy.o global.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./semilog.f90
set_flow_fractions.o: ./set_flow_fractions.f90 selective1usgs.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./set_flow_fractions.f90
shadec.o: ./shadec.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./shadec.f90
shading.o: ./shading.f90 shadec.o global.o gdayc.o surfhe.o geomc.o screenc.o logicc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./shading.f90
spillway_flow.o: ./spillway_flow.f90 structures.o global.o geomc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./spillway_flow.f90
spline.o: ./spline.f90 fishy.o global.o geomc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./spline.f90
structures.o: ./structures.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./structures.f90
surfhe.o: ./surfhe.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./surfhe.f90
tag124578.o: ./tag124578.f90 fishy.o global.o geomc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./tag124578.f90
tag369.o: ./tag369.f90 fishy.o global.o geomc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./tag369.f90
tdgas.o: ./tdgas.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./tdgas.f90
temperature.o: ./temperature.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o cemavars.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./temperature.f90
time_varying_data.o: ./time_varying_data.f90 global.o surfhe.o screenc.o tvdc.o logicc.o selwc.o structures.o namesc.o kinetic.o shadec.o main.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./time_varying_data.f90
total_dissolved_gas.o: ./total_dissolved_gas.f90 tdgas.o structures.o global.o main.o tvdc.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./total_dissolved_gas.f90
trans.o: ./trans.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./trans.f90
transport.o: ./transport.f90 global.o geomc.o tvdc.o trans.o logicc.o structures.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./transport.f90
tridiag.o: ./tridiag.f90 tridiag_v.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./tridiag.f90
tridiag_v.o: ./tridiag_v.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./tridiag_v.f90
tvdc.o: ./tvdc.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./tvdc.f90
twidth.o: ./twidth.f90 F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./twidth.f90
vgillnets.o: ./vgillnets.f90 global.o fishy.o geomc.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./vgillnets.f90
w2_dialog.o: ./w2_dialog.f90 
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./w2_dialog.f90
wallfunction.o: ./wallfunction.f90 global.o geomc.o eddy.o logicc.o main.o screenc.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./wallfunction.f90
waterbody.o: ./waterbody.f90 global.o geomc.o tvdc.o logicc.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./waterbody.f90
wetper.o: ./wetper.f90 F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./wetper.f90
whatjr.o: ./whatjr.f90 fishy.o global.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./whatjr.f90
#winmain.o: ./winmain.f90 F77KINDS.o
#	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./winmain.f90
withdrawal.o: ./withdrawal.f90 global.o geomc.o tvdc.o selwc.o logicc.o main.o screenc.o F77KINDS.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./withdrawal.f90
wqconstituents.o: ./wqconstituents.f90 main.o global.o namesc.o geomc.o logicc.o surfhe.o kinetic.o shadec.o eddy.o structures.o trans.o tvdc.o selwc.o gdayc.o screenc.o tdgas.o rstart.o macrophytec.o porosityc.o zooplanktonc.o tridiag_v.o cemavars.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./wqconstituents.f90
xsectional_area.o: ./xsectional_area.f90 global.o geomc.o main.o initialvelocity.o
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./xsectional_area.f90
zbrent1.o: ./zbrent1.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./zbrent1.f90
zooplanktonc.o: ./zooplanktonc.f90
	$(FC) $(FFLAGS) $(OTHERFLAGS) $(OTHER_FFLAGS) -c	./zooplanktonc.f90
SRC = ./twidth.f90 ./porosityc.f90 ./semilog.f90 ./cemabubblesturbulence.f90 ./hydroinout.f90 ./envirp.f90 ./deg_c.f90 ./fish.f90 ./cemacomputeturbidity.f90 ./cemabubblestransport.f90 ./selective1.f90 ./fishplot.f90 ./aerate.f90 ./endsimulation.f90 ./time_varying_data.f90 ./gas_transfer.f90 ./cemabubblesrelease.f90 ./output.f90 ./tag124578.f90 ./wallfunction.f90 ./selectiveinitusgs.f90 ./selectiveusgs.f90 ./tdgas.f90 ./fishy.f90 ./cemabubblesrisevelocity.f90 ./interflowf.f90 ./zbrent1.f90 ./total_dissolved_gas.f90 ./spline.f90 ./heat_exchange.f90 ./tridiag.f90 ./init.f90 ./outputa.f90 ./findnewbr.f90 ./density.f90 ./spillway_flow.f90 ./ludcmp.f90 ./kinetics.f90 ./waterbody.f90 ./surfhe.f90 ./part_transport.f90 ./shadec.f90 ./fishhabitat.f90 ./lubksb.f90 ./rtbis.f90 ./selective.f90 ./balances.f90 ./input.f90 ./initgeom.f90 ./calcfric.f90 ./gasbubblesformation.f90 ./geomc.f90 ./cdfunc.f90 ./cemacalculaterisevelocity.f90 ./set_flow_fractions.f90 ./grid_area2.f90 ./xsectional_area.f90 ./findbranch.f90 ./selwc.f90 ./wqconstituents.f90 ./bioenergetics.f90 ./gasdev.f90 ./eddy.f90 ./screen_update.f90 ./restart_output.f90 ./calculate_az.f90 ./read_fish_data.f90 ./pipe_flow_initialize.f90 ./gdayc.f90 ./fimpbr.f90 ./porosity.f90 ./cemafftlayercode.f90 ./w2_dialog.f90 ./transport.f90 ./shading.f90 ./interconst.f90 ./gate_flow.f90 ./cema_w2_input.f90 ./barea.f90 ./temperature.f90 ./grid_area1.f90 ./acoustics.f90 ./structures.f90 ./main.f90 ./mannings_eqn.f90 ./gregorian_date.f90 ./normal_depth.f90 ./selective1usgs.f90 ./trans.f90 ./macrophytec.f90 ./gridplot.f90 ./tag369.f90 ./whatjr.f90 ./withdrawal.f90 ./envirpmod.f90 ./initial_u_velocity.f90 ./cemavars.f90 ./initializebedconsolidationfiles.f90 ./initcond.f90 ./zooplanktonc.f90 ./wetper.f90 ./cemasedimentmodelw2.f90 ./screenc.f90 ./vgillnets.f90 ./depthcrit.f90 ./initialvelocity.f90 ./cemadisgasphasedistribution.f90 ./lateral_velocity.f90 ./layeraddsub.f90 ./cemabubblesreleaseturbulence.f90 ./readbedconsolidationfiles.f90 ./cemabubbwattransfer.f90 ./fishoutput.f90 ./deg_f.f90 ./selectiveinit.f90 ./cemaoutputroutines.f90 ./lin_sys.f90 ./initial_water_level.f90 ./enable.f90 ./tridiag_v.f90 ./F77KINDS.F90 ./cemasedimentdiagenesis.f90 ./tvdc.f90 ./deallocate_cema.f90 ./random.f90 ./open_channel_initialize.f90 ./macrophyte_friction.f90 ./global.f90 ./namesc.f90 ./logicc.f90 ./run_w2.f90 ./rstart.f90 ./kinetic.f90
OBJ = twidth.o porosityc.o semilog.o cemabubblesturbulence.o hydroinout.o envirp.o deg_c.o fish.o cemacomputeturbidity.o cemabubblestransport.o selective1.o fishplot.o aerate.o endsimulation.o time_varying_data.o gas_transfer.o cemabubblesrelease.o output.o tag124578.o wallfunction.o selectiveinitusgs.o selectiveusgs.o tdgas.o fishy.o cemabubblesrisevelocity.o interflowf.o zbrent1.o total_dissolved_gas.o spline.o heat_exchange.o tridiag.o init.o outputa.o findnewbr.o density.o spillway_flow.o ludcmp.o kinetics.o waterbody.o surfhe.o part_transport.o shadec.o fishhabitat.o lubksb.o rtbis.o selective.o balances.o input.o initgeom.o calcfric.o gasbubblesformation.o geomc.o cdfunc.o cemacalculaterisevelocity.o set_flow_fractions.o grid_area2.o xsectional_area.o findbranch.o selwc.o wqconstituents.o bioenergetics.o gasdev.o eddy.o screen_update.o restart_output.o calculate_az.o read_fish_data.o pipe_flow_initialize.o gdayc.o fimpbr.o porosity.o cemafftlayercode.o w2_dialog.o transport.o shading.o interconst.o gate_flow.o cema_w2_input.o barea.o temperature.o grid_area1.o acoustics.o structures.o main.o mannings_eqn.o gregorian_date.o normal_depth.o selective1usgs.o trans.o macrophytec.o gridplot.o tag369.o whatjr.o withdrawal.o envirpmod.o initial_u_velocity.o cemavars.o initializebedconsolidationfiles.o initcond.o zooplanktonc.o wetper.o cemasedimentmodelw2.o screenc.o vgillnets.o depthcrit.o initialvelocity.o cemadisgasphasedistribution.o lateral_velocity.o layeraddsub.o cemabubblesreleaseturbulence.o readbedconsolidationfiles.o cemabubbwattransfer.o fishoutput.o deg_f.o selectiveinit.o cemaoutputroutines.o lin_sys.o initial_water_level.o enable.o tridiag_v.o F77KINDS.o cemasedimentdiagenesis.o tvdc.o deallocate_cema.o random.o open_channel_initialize.o macrophyte_friction.o global.o namesc.o logicc.o run_w2.o rstart.o kinetic.o
clean: neat
	-rm -f .a.out.cppdefs $(OBJ) a.out

neat:
	-rm -f $(TMPFILES)

TAGS: $(SRC)
	etags $(SRC)

tags: $(SRC)
	ctags $(SRC)

a.out: $(OBJ)
	$(LD) $(OBJ) -o a.out  $(LDFLAGS)

