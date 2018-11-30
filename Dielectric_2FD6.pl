#!perl

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);
use List::Util qw[min max];
use Math::Trig;


#------------------------------------------------------------------------------------------------------------------------------------------
#                                               Carregando arquivo PDB
#------------------------------------------------------------------------------------------------------------------------------------------

my $doc = $Documents{"2FD6-uPA-uPAR_addH_opt.xsd"};  #Fornecer o nome do arquivo .xsd ("nome_do_arquivo.xsd")

#------------------------------------------------------------------------------------------------------------------------------------------
#                                             Fazendo copia do PDB sem H2O
#------------------------------------------------------------------------------------------------------------------------------------------

my $doc1=no_water($doc);

#------------------------------------------------------------------------------------------------------------------------------------------
#                                              Criando tabela Dielectric.std
#------------------------------------------------------------------------------------------------------------------------------------------

my $table = Documents->New("Dielectric_A1_CTLA-4.std");

$table -> ColumnHeading(0) = "Residue_uPA"; 
$table -> ColumnHeading(1) = "Residue_uPAR"; 
$table -> ColumnHeading(2) = "Dielectric constant"; 
$table -> ColumnHeading(3) = "Distance"; 
$table -> ColumnHeading(4) = "# de pontos(%)"; 
$table -> ColumnHeading(5) = "Xc_sphere"; 
$table -> ColumnHeading(6) = "Yc_sphere"; 
$table -> ColumnHeading(7) = "Zc_sphere"; 
$table -> ColumnHeading(8) = "R_sphere"; 
$table -> ColumnHeading(9) = "Error(%)"; 
$table -> UpdateViews;

#------------------------------------------------------------------------------------------------------------------------------------------
#                                               Criando listas A1,A2 e B
#------------------------------------------------------------------------------------------------------------------------------------------

my @list_A;
my @list_B;

my $new_chain;
my $id_residue=-1;
foreach my $residue (@{$doc1->SubUnits}){
	$id_residue=$id_residue+1;
	
	my $chain=$residue->Ancestors->Chain->Name;
	my $name=$residue->Name;
	
	if($chain eq "A"){push(@list_A,$id_residue)}
	
	if($chain eq "U"){push(@list_B,$id_residue)}
}

#------------------------------------------------------------------------------------------------------------------------------------------
#                                             Parametros da rede cubica
#------------------------------------------------------------------------------------------------------------------------------------------

pdb_quadrante1($doc1);          #Transladando proteína para o primeiro quadrante

my $r_cube=3;                   #Definindo largura dos cubos pequenos
my $dmax=length_cube($doc1);
my $L=int($dmax/$r_cube)+1;     #Definindo a quantidade de cubos por "linha"

#------------------------------------------------------------------------------------------------------------------------------------------
#                                             Construindo rede cubica
#------------------------------------------------------------------------------------------------------------------------------------------

my @cube_array=cube_space($doc1,$r_cube,$L);   #L -> largura da rede cubica

#------------------------------------------------------------------------------------------------------------------------------------------
#                                  Calculando a constante dieletrica para interacoes A1 e B
#------------------------------------------------------------------------------------------------------------------------------------------

my $count=0;
my $row=-1;
foreach my $id1 (@list_A){	
	
	my $res1=$doc1->SubUnits->Item($id1)->Name;
	
	foreach my $id2 (@list_B){	
	
	my $res2=$doc1->SubUnits->Item($id2)->Name;
	my($dist,$atom1,$atom2)=res_dist($doc1,$id1,$id2);
	
		if($dist<=10){
			$row=$row+1;
			$count=$count+1;			
			
			my ($x1,$y1,$z1)=centroid($doc1,$id1);
			my ($x2,$y2,$z2)=centroid($doc1,$id2);
			my ($x_sphere,$y_sphere,$z_sphere,$R_sphere)=sphere($doc1,$id1,$x1,$y1,$z1,$id2,$x2,$y2,$z2);
			
		#-------------------------------------------------------------------------
		#      Definindo variaveis associadas ao cálculo da media da funcao
		#-------------------------------------------------------------------------
			
			my $fm=0;
			my $fm2=0;
			
			my $fm_average;
			my $fm2_average;
			
			my $function;
			my $sum_function=0;
			
			my $nbox=0;
			my $box_size=1000;
			
			my $error=10000;
			my $npoints=0;
				
		#-------------------------------------------------------------------------
		#      						Grid aleatorio
		#-------------------------------------------------------------------------
		
			while($npoints<7000 or $error>2.5){			
				
				$npoints=$npoints+1;
				
				my $pi=3.14159265;
				
				my $theta=acos(1-2*rand());
				my $phi=2*$pi*rand();
				my $r=$R_sphere*rand();
			
				my $x= $x_sphere + $r * sin($theta) * cos($phi);
				my $y= $y_sphere + $r * sin($theta) * sin($phi);
				my $z= $z_sphere + $r * cos($theta);
			
				my $id_cube=id_cube($r_cube,$L,$x,$y,$z);	
				
				$function=dielectric_function_app($doc1,$id_cube,$L,$x,$y,$z,\@cube_array);
				$sum_function=$sum_function+$function;
				
				#printf "%10d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f \n",$npoints,$x,$y,$z,$r,$R_sphere,$function,$error;					
			
				if($npoints % $box_size ==0){
					$nbox=$nbox+1;
					
					my $box_average = $sum_function/$box_size;
					
					$fm  = $fm  + $box_average;	
					$fm2 = $fm2 + $box_average**2;
		
					$fm_average=$fm/$nbox;
					$fm2_average=$fm2/$nbox;
					
					if($nbox>1){
						$error=sqrt( ( $fm2_average - ($fm_average**2) ) /$nbox);
						$error=100*1.96*$error/$fm_average;	
					}
				
					#printf "%10d %10.4f %10.4f %10.4f \n",$npoints,$box_average,$fm_average,$error;	
				
					$sum_function=0;
				}
					
			}
	
			printf "%10s %10s %10.4f %10.4f %10d %10.4f %10.4f %10.4f %10.4f %10.4f \n",$res1,$res2,$fm_average,$dist,$npoints,$x_sphere,$y_sphere,$z_sphere,$R_sphere,$error;	
				
			
			$table->cell($row,0)=$res1;
			$table->cell($row,1)=$res2;
			$table->cell($row,2)=$fm_average;
			$table->cell($row,3)=$dist;
			$table->cell($row,4)=$npoints;
			$table->cell($row,5)=$x_sphere;
			$table->cell($row,6)=$y_sphere;
			$table->cell($row,7)=$z_sphere;	
			$table->cell($row,8)=$R_sphere;
			$table->cell($row,9)=$error;
			$table->UpdateViews;
		
		}
	
	}
}


#-------------------------------------------------------------------------
#                           sub no_water
#-------------------------------------------------------------------------

sub no_water{
	my $doc_aux=@_[0];

	my $doc1=Documents->New("PDB_no_water.xsd");
	$doc1->CopyFrom($doc_aux);
	
	my $water_chain;
	
	foreach my $chain (@{$doc1->Chains}){
		my $residue=$chain->SubUnits->Item(0)->Name;
		if(substr($residue,0,3) eq "HOH"){
			$water_chain=$chain->Name;
			$doc1->Chains->DisplayRange($water_chain)->Delete;
		}
	}

	return $doc1;
}


#-------------------------------------------------------------------------
#                           sub res_caps
#-------------------------------------------------------------------------

sub res_caps{

	my $doc_aux=@_[0];
	my $id_residue=@_[1];
	my $lim_down=@_[2];
	my $lim_up=@_[3];
	my $caps_number=@_[4];
	
	my $i;
	my @caps;
	
	
	my $ncaps=-1;
	for ($i=1;$i<=$caps_number;$i++){
		
#            -------------------------
#                   Caps à direita
#            -------------------------
		
		my $id_cap_r = $id_residue + $i;
		
		if($id_cap_r <= $lim_up){
			$ncaps=$ncaps+1;
			$caps[$ncaps] = $id_cap_r;
		}
	
#            -------------------------
#                   Caps à esquerda
#            -------------------------
		
		my $id_cap_l = $id_residue - $i;

		if($id_cap_l >= $lim_down){
			$ncaps=$ncaps+1;
			$caps[$ncaps] = $id_cap_l;
		}
	
	}
	
	return ($ncaps,@caps);

}

#-------------------------------------------------------------------------
#                           sub res_dist
#-------------------------------------------------------------------------

sub res_dist{

	my $doc_aux=@_[0];
	my $i=@_[1];
	my $j=@_[2];
	
	my $dist_min=1000;
	my $atom1_name;
	my $atom2_name;
	
	foreach my $atom1 (@{$doc_aux->SubUnits->Item($i)->Atoms}){
			foreach my $atom2 (@{$doc_aux->SubUnits->Item($j)->Atoms}){
				my $dist = sqrt ( ($atom1->X - $atom2->X)**2 + ($atom1->Y - $atom2->Y)**2 + ($atom1->Z - $atom2->Z)**2 );
				if ($dist<$dist_min){
						$dist_min=$dist;
						$atom1_name=$atom1->Name;
						$atom2_name=$atom2->Name;
				}	
			}
		}

	return ($dist_min,$atom1_name,$atom2_name);

}

#-------------------------------------------------------------------------
#                       sub pdb_quadrante1
#-------------------------------------------------------------------------

sub pdb_quadrante1{
	
	my $doc_aux=@_[0];

	my $xmin=0;
	my $ymin=0;
	my $zmin=0;
	
	my @atom_set;
	foreach my $atom (@{$doc_aux->Atoms}){
		push(@atom_set,$atom);
		
		my $x=$atom->X;
		my $y=$atom->Y;
		my $z=$atom->Z;
			
		if($x<$xmin){$xmin=$x}
		if($y<$ymin){$ymin=$y}
		if($z<$zmin){$zmin=$z}

	}

	my $dx=abs($xmin);
	my $dy=abs($ymin);
	my $dz=abs($zmin);


	$doc_aux->CreateSet("all_atoms", \@atom_set); #Criando conjunto com todos os atomos 
	$doc_aux->Sets("all_atoms")->Translate(Point(X => $dx, Y => $dy, Z => $dz));

}

#-------------------------------------------------------------------------
#                          sub length_cube
#-------------------------------------------------------------------------

sub length_cube{
	my $doc_aux=@_[0];
	
	my $Lmax=-10000;
	foreach my $atom (@{$doc_aux->Atoms}){
		my $x=$atom->X;
		my $y=$atom->Y;
		my $z=$atom->Z;
	
		if(abs($x)>$Lmax){$Lmax=abs($x)}
		if(abs($y)>$Lmax){$Lmax=abs($y)}
		if(abs($z)>$Lmax){$Lmax=abs($z)}
	
	}

	return $Lmax;	
}

#-------------------------------------------------------------------------
#                      sub dielectric_function_app
#-------------------------------------------------------------------------

sub dielectric_function_app{
	
	my $doc_aux=@_[0];
	my $id=@_[1];
	my $L=@_[2];
	
	my $x=@_[3];
	my $y=@_[4];
	my $z=@_[5];
	
	my @cube_array=@{$_[6]};
	
	my @v;

	my $N=$L*$L;

	$v[0]=$id-1;
	$v[1]=$id+1;
	$v[2]=$id+$L;
	$v[3]=$id-$L;
	$v[4]=$id-$L-1;
	$v[5]=$id-$L+1;
	$v[6]=$id+$L-1;
	$v[7]=$id+$L+1;
	
	
	$v[8]=$id+$N;
	$v[9]=$id+$N+1;
	$v[10]=$id+$N-1;
	$v[11]=$id+$N+$L;
	$v[12]=$id+$N-$L;
	$v[13]=$id+$N+$L+1;
	$v[14]=$id+$N+$L-1;
	$v[15]=$id+$N-$L+1;
	$v[16]=$id+$N-$L-1;
	

	$v[17]=$id-$N;
	$v[18]=$id-$N+1;
	$v[19]=$id-$N-1;
	$v[20]=$id-$N+$L;
	$v[21]=$id-$N-$L;
	$v[22]=$id-$N+$L+1;
	$v[23]=$id-$N+$L-1;
	$v[24]=$id-$N-$L+1;
	$v[25]=$id-$N-$L-1;
	
	$v[26]=$id;
	
	my $sigma=0.93;           #parametro da distribuição gaussiana considerada por Li et al.
	my $epslon_in=4;          #valor de referencia da constante dieletrica no centro do átomo
	my $epslon_out=80;        #valor de referencia da constante dieletrica no solvente (agua)

	
	my $density_mol=1;
	foreach my $id_aux (@v){	
		
		if($id_aux >=1 and $id_aux <= ($L*$L*$L)){
			my $aux1=$#{$cube_array[$id_aux]};
						
			for(my $i = 0; $i <= $#{$cube_array[$id_aux]}; $i++){
				my $aux2=$cube_array[$id_aux][$i];
									
				my $atom=$doc_aux->Atoms->Item($aux2);
				my $VDWR=$atom->VDWRadius;
				my $dist=sqrt(  ($x - $atom->X)**2   +  ($y - $atom->Y)**2  +  ($z - $atom->Z)**2  );
					
				my $density_atom=exp(  -$dist**2/ ( ($sigma**2)*($VDWR**2) )  );			
				$density_mol=$density_mol*(1-$density_atom);                                               		
								
			}
				
		}
	
	}


	$density_mol=1-$density_mol;
	my $dielectric   = $density_mol*$epslon_in + (1-$density_mol)*$epslon_out;
	
	return $dielectric;
}

#-------------------------------------------------------------------------
#                           sub sphere
#-------------------------------------------------------------------------
# O raio da esfera e definido como a maior distancia residuo/ligante para o centroide residuo-ligante		

sub sphere{

	my $doc_aux=@_[0];	
	
	my $id1=@_[1];
	my $x1=@_[2];
	my $y1=@_[3];
	my $z1=@_[4];
	
	my $id2=@_[5];
	my $x2=@_[6];
	my $y2=@_[7];
	my $z2=@_[8];
		
	my $x_sphere;
	my $y_sphere;
	my $z_sphere;
	my $R_sphere;
	
	$x_sphere=($x1+$x2)/2;	
	$y_sphere=($y1+$y2)/2;
	$z_sphere=($z1+$z2)/2;	
		
	my $max_dist=0;
	foreach my $atom (@{$doc_aux->SubUnits->Item($id1)->Atoms}){
		$R_sphere = sqrt(($atom->X-$x_sphere)**2 + ($atom->Y-$y_sphere)**2 + ($atom->Z-$z_sphere)**2);	
		if($R_sphere>$max_dist){$max_dist=$R_sphere;}
	}	
	
	foreach my $atom (@{$doc_aux->SubUnits->Item($id2)->Atoms}){
		$R_sphere = sqrt(($atom->X-$x_sphere)**2 + ($atom->Y-$y_sphere)**2 + ($atom->Z-$z_sphere)**2);	
		if($R_sphere>$max_dist){$max_dist=$R_sphere;}				
	}
	
	$R_sphere=$max_dist;

	return ($x_sphere,$y_sphere,$z_sphere,$R_sphere);
}

#-------------------------------------------------------------------------
#                           sub centroid
#-------------------------------------------------------------------------

sub centroid {

	my $doc_aux=@_[0];
	my $id_residue=@_[1];
	
	my $centroid=$doc_aux->CreateCentroid($doc_aux->SubUnits->Item($id_residue)->Atoms);    #Calculando centro geometrico do ligante
	$centroid->IsWeighted = "No";    						                                #Opcao por centro geometrico ao inves do centro de massa (media ponderada pelas massas atomicas)
	
	my $center=$centroid->CentroidXYZ;
	my $x=$center->X;
	my $y=$center->Y;
	my $z=$center->Z;

	return ($x,$y,$z);

}

#-------------------------------------------------------------------------
#                        sub id_cube
#-------------------------------------------------------------------------

sub id_cube{

	my $r=@_[0];
	my $L=@_[1];
	my $x=@_[2];
	my $y=@_[3];
	my $z=@_[4];
	
	my $idx=int($x/$r)+1;
	my $idy=int($y/$r)+1;
	my $idz=int($z/$r)+1;

	my $id = $idx + $L*($idy-1) + $L*$L*($idz-1);
	
	return $id;

}

#-------------------------------------------------------------------------
#                       sub cube_space
#-------------------------------------------------------------------------

sub cube_space{

	my $doc_aux=@_[0];
	my $r=@_[1];
	my $L=@_[2];
	
	my $id;
	my $idx;
	my $idy;
	my $idz;
	
	my @cube_array;
	
	my $id_atom=-1;
	foreach my $atom (@{$doc_aux->Atoms}){
		
		$id_atom=$id_atom+1;
		
		my $x=$atom->X;
		my $y=$atom->Y;
		my $z=$atom->Z;
				
		$idx=int($x/$r)+1;
		$idy=int($y/$r)+1;
		$idz=int($z/$r)+1;
		
		$id = $idx + $L*($idy-1) + $L*$L*($idz-1);
		
		push(@{$cube_array[$id]},$id_atom);
		
		#printf "%10d %10d %10d %10d %10.2f %10.2f %10.2f %10d %10d %10d \n",$r,$L,$id,$id_atom,$x,$y,$z,$idx,$idy,$idz;
		
	}

	return @cube_array;
}





