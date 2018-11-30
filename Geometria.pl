################Otimização dos hidrogênios e cálculo da energia da molécula###########################
#Descrição: Esse script executa a seguinte tarefa:
#1 - Otimiza a geometria dos hidrogênios
##################################Atulização:30-08-2018###############################################

#!perl
use strict;
use Getopt::Long;
use MaterialsScript qw(:all);
use Cwd;


my $dir="C:/Users/Pesquisa/Google Drive/uPA-uPAR_Files/Documents";


opendir(diretorio, "$dir");
my @lista = grep { /\.xsd$/ } readdir(diretorio);
@lista = grep(/^2FD6/, @lista);                    
closedir(diretorio);


foreach my $key (@lista) {
	print "processing file key $key\n";
	
	my $doc = $Documents{$key};
	
	Tools->BondCalculation->ChangeSettings(Settings(ResonantBondRepresentation => "Resonant",));
	$doc -> CalculateBonds; 
	
	foreach my $atom (@{$doc->Atoms}) {
		if ($atom->ElementSymbol ne "H") { $atom->Fix("XYZ");}
	}
	

	Modules->Forcite->ChangeSettings([CurrentForcefield => "Universal"]); 
	Modules->Forcite->ChangeSettings([Quality => "Ultra-fine"]); 
	Modules->Forcite->ChangeSettings([ChargeAssignment => "Forcefield assigned"]); 
	Modules->Forcite->ChangeSettings([NonPeriodicElectrostaticSummationMethod => "Atom based"]);
	Modules->Forcite->ChangeSettings([NonPeriodicvdWSummationMethod => "Atom based"]); 
	Modules->Forcite->ChangeSettings([NonPeriodicHBondSummationMethod => "Atom based"]); 
	Modules->Forcite->ChangeSettings([MaxIterations => "50000"]); 
	Modules->Forcite->ChangeSettings([MaxChargeIterations => "50000"]); 
	Modules->Forcite->ChangeSettings([OptimizationAlgorithm => "Smart"]); 
	Modules->Forcite->ChangeSettings([NonPeriodicvdWAtomTruncationMethod => "None"]);  
	
	Modules -> Forcite -> GeometryOptimization -> run($doc); 

	foreach my $atom (@{$doc->Atoms}) {
		if ($atom->ElementSymbol ne "H") { $atom->Unfix("XYZ");}
	}

}



