#!perl

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);
use List::Util qw[min max];
use Cwd;


my $table1 = Documents->New("Energy.std");
my $table2 = $Documents{"Residues_2FD6_uPA_uPAR.std"};                #Carregando tabela(existente) "Residues.std"


$table1 -> ColumnHeading(0)="R_uPA"; 
$table1 -> ColumnHeading(1)="R_uPAR";

$table1 -> ColumnHeading(2)="Energy A (kcal/mol)"; 
$table1 -> ColumnHeading(3)="Energy_B (kcal/mol)"; 
$table1 -> ColumnHeading(4)="Energy_C (kcal/mol)"; 
$table1 -> ColumnHeading(5)="Energy_D (kcal/mol)"; 

$table1 -> ColumnHeading(6)="Distance (Angstrons)";
$table1 -> ColumnHeading(7)="Dielectric Constant";
$table1 -> ColumnHeading(8)="Interaction Energy (kcal/mol)"; 

$table1 -> ColumnHeading(9)="R_uPA Type";
$table1 -> ColumnHeading(10)="R_uPA Charge";
$table1 -> ColumnHeading(11)="R_uPA Chain";
$table1 -> ColumnHeading(12)="R_uPA Atom";

$table1 -> ColumnHeading(13)="R_uPAR Type";
$table1 -> ColumnHeading(14)="R_uPAR Charge";
$table1 -> ColumnHeading(15)="R_uPAR Atom";
$table1 -> ColumnHeading(16)="R_uPAR Chain";

$table1 -> ColumnHeading(17)="Interaction Energy (Ha)";
$table1 -> ColumnHeading(18)="File Name A";
$table1 -> ColumnHeading(19)="File Name B";
$table1 -> ColumnHeading(20)="File Name C";
$table1 -> ColumnHeading(21)="File Name D";
$table1->UpdateViews;

my $E_A;
my $E_B;
my $E_C;
my $E_D;

my $energy=0;

my $dielectric_function="yes";
my $constante;
my $table3;

if($dielectric_function eq "no"){
	$constante=40.00;
}
else{
	$table3 = $Documents{"Dielectric_2FD6.std"};
}

#my $cwd="Perl%20Script_Files\\Documents";
my $cwd="C:/Users/Pesquisa/Google Drive/uPA-uPAR_Files/Documents/2FD6/MFCC_2FD6/MFCC_2FD6 Script";

opendir my $folder, $cwd or die "Couldn't open folder";
my @allfiles = grep {!/^\.\.?$/} readdir $folder;                #Armazenando o nome de todos os arquivos
@allfiles = grep(/^U1_D/, @allfiles);							 #Selecionando só os nomes que comecam com "Dist_"	

my $row1=-1;
foreach my $file (@allfiles) {
 
   my $item=substr($file, 0, index($file, '.'));    #extrai todos os caracteres anteriores ao ponto (retira a extensão .xsd)
   $item=substr($item,-1,1);                        # Identifica se é arquivo A, B, C ou D 
   if($item eq "A"){$row1=$row1+1;}
   
#------------------------------------------------------------------------
#		    Escrevendo nomes dos arquivos (A,B,C,D,E e F) 
#------------------------------------------------------------------------

   if($item eq "A"){$table1->cell($row1,18)=$file;}
   if($item eq "B"){$table1->cell($row1,19)=$file;}
   if($item eq "C"){$table1->cell($row1,20)=$file;}
   if($item eq "D"){$table1->cell($row1,21)=$file;}
	  
#------------------------------------------------------------------------
#	   Dividindo nomes dos arquivos em campos delimitados por "_"
#------------------------------------------------------------------------

   my @words = split /_/,$file;
   
   my $residue1=$words[4];
   my $residue2=$words[5];
  
   my $residue1_type=substr($residue1,0,3);
   my $residue2_type=substr($residue2,0,3);
   
   my $line; 

#------------------------------------------------------------------------
#				lendo tabela do MFCC (Residues.std)
#------------------------------------------------------------------------
	my $i;
	my $nline=115;
	my $dist;
	my $chain1;
	my $chain2;
	my $atom1;
	my $atom2;
	my $charge1;
	my $charge2;
	for ($i=0;$i<$nline;$i++){   #i < numero de linhas da tabela
		if($residue1 eq $table2->cell($i,0) and $residue2 eq $table2->cell($i,1)){		
			$dist=$table2->cell($i,2);
			$chain1=$table2->cell($i,3);	
			$atom1=$table2->cell($i,4);
			$atom2=$table2->cell($i,5);
			$chain2=$table2->cell($i,6);
			$charge1=$table2->cell($i,7);
			$charge2=$table2->cell($i,8);
			last;
		}
	}

  
	$table1->cell($row1,0)=$residue1;
	$table1->cell($row1,1)=$residue2;
	 
	$table1->cell($row1,6)=$dist;
	$table1->cell($row1,7)=$constante;
	  
	$table1->cell($row1,9)=$residue1_type;
	$table1->cell($row1,10)=$charge1;
	$table1->cell($row1,11)=$chain1;
	$table1->cell($row1,12)=$atom1;
		   
	$table1->cell($row1,13)=$residue2_type;
	$table1->cell($row1,14)=$charge2;
	$table1->cell($row1,15)=$atom2;
	$table1->cell($row1,16)=$chain2;	   

#------------------------------------------------------------------------
#	               Atribuindo constante dieletrica
#------------------------------------------------------------------------	 

   $nline=718;
   if($dielectric_function eq "yes"){
	   for ($i=0;$i<$nline;$i++){   #i < numero de linhas da tabela
			if($residue1 eq $table3->cell($i,0) and $residue2 eq $table3->cell($i,1)){		
				$constante=$table3->cell($i,2);
				last;
			}
		}
	}

#---------------------------------------------------------------------
#                     Carregando Arquivo
#---------------------------------------------------------------------
	
	my $doc = $Documents{"$file"};       #Carregando o arquivo .xsd

#---------------------------------------------------------------------
#              Calculo da carga formal (Formal charge)
#---------------------------------------------------------------------

	my $formalChargeSum = 0;
	foreach my $atom (@{$doc->Atoms}){
		my $formalCharge = $atom->FormalCharge->Value;
		#printf "Formal charge on atom %s = %f/%f\n", $atom->Name, $formalCharge;
		$formalChargeSum += $formalCharge;
	}

#---------------------------------------------------------------------
#                Configurando Dmol  (Set up Dmol)
#---------------------------------------------------------------------

	my $dmol3 = Modules->DMol3;
			
		        $dmol3->ChangeSettings([Quality => "Fine",
		        			
		        					UseSymmetry => "No",
			                        
			                		ElectronicQuality => "Fine",
			                        
			                        Charge => "$formalChargeSum",
			                        
			                        Multiplicity => "Auto",
			                        
			                        SpinUnrestricted => "Yes",
			                        
			                        TheoryLevel => "GGA",		#LDA ou GGA.
						
									# Usar LocalFunciontal p/ LDA
									# Usar NonLocalFunctional p/ GGA
						
			                        #LocalFunctional => "PWC",	#Comentar caso use GGA.
	
									NonLocalFunctional=> "PBE",	#Apenas para GGA.
						
									UseDFTD => "Yes",
						
									DFTDMethod => "TS",	    #OBS, TS ou Grimme.
	
			                        Basis => "DNP+", 	        #MIN, DN, DND, DNP, MIXED.
			                        
			                        BasisFile => "4.4",	        #"3.5" ou "4.4"
			                        
			                        MaximumSCFCycles => "1000",
			                        
			                        UseSmearing => "YES",
			                        
			                        Smearing => "0.005",
			                      
			                        UseDIIS => "YES",
			                        
			                        CutoffType => "Custom",
			       					
			       					AtomCutoff => "5.5",
			                        
			                        CoreTreatment => "All Electron", 	# "All Electron"
			                        									# "Effectove core potentials"
			                        									# "All electrons relativistic"
			                        									# "DFT semi-core pseudoposts"
			                        UseCosmo => "YES", 					#Solvente
			                        
			                        #CosmoSolvent => "Water",		# "Water", "Dimethyl sulfoxide"
			                        
			                        SolventDielectric => "$constante",	# 1 a 10000 (Solvent dielectric constat)
			                        
			                                  
									#PopulationAnalysis => "Yes",
	]);


#---------------------------------------------------------------------
#                     Executando Dmol  (Run Dmol)
#---------------------------------------------------------------------
	
	my $output=$dmol3->Energy->Run($doc);
	my $FinalEnergy=$output->TotalEnergy;
 	 
 	#my $FinalEnergy=20;
 	 	
#---------------------------------------------------------------------
#           Escrevendo as energias na tabela do Dmol (Energy.std)
#---------------------------------------------------------------------

	if ($item eq "A"){
		$table1->cell($row1,2)=$FinalEnergy;
		$E_A=$FinalEnergy;
		$table1->UpdateViews;
	}
	
	if ($item eq "B"){
		$table1->cell($row1,3)=$FinalEnergy;
		$E_B=$FinalEnergy;
		$table1->UpdateViews;
	}
	
	if ($item eq "C"){
		$table1->cell($row1,4)=$FinalEnergy;
		$E_C=$FinalEnergy;
		$table1->UpdateViews;
	}
	
	if ($item eq "D"){
		$table1->cell($row1,5)=$FinalEnergy;
		$E_D=$FinalEnergy;
		$energy= $E_A-$E_B-$E_C+$E_D;
		$table1->cell($row1,8)=$energy;
		$table1->cell($row1,17)=$energy/627.509;	
		$table1->UpdateViews;
		printf "%20s %10d %10d %20s %20.4f %20.2f %20.5f %20.5f %20.5f %20.5f %20.5f \n",$residue1,$charge1,$charge2,$residue2,$dist,$constante,$E_A,$E_B,$E_C,$E_D,$energy;
	}
	
	$doc->Close;
}

