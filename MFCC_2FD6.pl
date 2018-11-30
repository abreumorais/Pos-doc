#!perl

use strict;
use Getopt::Long;
use MaterialsScript qw(:all);

#------------------------------------------------------------------------------------------------------------------------------------------
#                                                  Carregando PDB 
#------------------------------------------------------------------------------------------------------------------------------------------

my $doc = $Documents{"2FD6-uPA-uPAR_addH_opt.xsd"};  #Fornecer o nome do arquivo .xsd ("nome_do_arquivo.xsd")


#------------------------------------------------------------------------------------------------------------------------------------------
#                                             Criando tabela Residues.std
#------------------------------------------------------------------------------------------------------------------------------------------

my $table2 = Documents->New("Residues_2FD6_uPA_uPAR.std");
$table2 -> ColumnHeading(0) = "R_uPA"; 
$table2 -> ColumnHeading(1) = "R_uPAR"; 
$table2 -> ColumnHeading(2) = "Distance"; 
$table2 -> ColumnHeading(3) = "R_uPA chain"; 
$table2 -> ColumnHeading(4) = "R_uPA atom"; 
$table2 -> ColumnHeading(5) = "R_uPAR atom"; 
$table2 -> ColumnHeading(6) = "R_uPAR chain"; 
$table2 -> ColumnHeading(7) = "R_uPA  charge";
$table2 -> ColumnHeading(8) = "R_uPAR charge";

#------------------------------------------------------------------------------------------------------------------------------------------
#                           Alterando nome dos residuos e transladando PDB para primeiro quadrante
#------------------------------------------------------------------------------------------------------------------------------------------

$doc=change_res_name($doc);
$doc=pdb_quadrante1($doc);
my $doc1=no_water($doc);

#------------------------------------------------------------------------------------------------------------------------------------------
#                                              Parametros da rede cubica
#------------------------------------------------------------------------------------------------------------------------------------------

my $dmax=length_cube($doc);
my $r_cube=12;                   #largura dos cubos pequenos
my $L=int($dmax/$r_cube)+1;      #Quantidade de cubos por "linha"

#------------------------------------------------------------------------------------------------------------------------------------------
#                                  Construindo rede cubica e mapeando moleculas de agua
#------------------------------------------------------------------------------------------------------------------------------------------

my @cube_array=cube_space($doc,$r_cube,$L);       

#------------------------------------------------------------------------------------------------------------------------------------------
#                                         Associando moléculas de água aos residuos
#------------------------------------------------------------------------------------------------------------------------------------------

my @water_array=water_array($doc1,$doc,$r_cube,$L,\@cube_array);

#------------------------------------------------------------------------------------------------------------------------------------------
#                                             Identificando cadeias e residuos
#------------------------------------------------------------------------------------------------------------------------------------------

my $id=-1;
foreach my $chain (@{$doc1->Chains}){
	$id=$id+1;
	my $name=$chain->Name;
	my $nres=$chain->SubUnits->Count;
	print "$id,$name,$nres\n";
}	

#------------------------------------------------------------------------------------------------------------------------------------------
#                                                     Criando listas A e B
#------------------------------------------------------------------------------------------------------------------------------------------

my @list_A;
my @list_U;

my $countA=0;
my $countU=0;

my $lim1_A;
my $lim2_A;

my $lim1_U;
my $lim2_U;
	
my $id=-1;
foreach my $residue (@{$doc1->SubUnits}){
	$id=$id+1;
	my $chain=$residue->Ancestors->Chain->Name;
	my $name=$residue->Name;
	
	my $charge1=formal_charge($doc1,$id);
	print "$id,$name,$charge1 \n";
	
	if($chain eq "A"){
		$countA=$countA+1;
		if($countA==1){$lim1_A=$id}
		$lim2_A=$id;
		push(@list_A,$id);
	}
	
	if($chain eq "U"){
		
		$countU=$countU+1;
		if($countU==1){$lim1_U=$id}
		$lim2_U=$id;
		push(@list_U,$id);
	}
				
}

print "$lim1_A,$lim2_A,$lim1_U,$lim2_U \n";

print "@list_A \n";
print "@list_U \n";

#------------------------------------------------------------------------------------------------------------------------------------------
#                                                          MFCC U1-A
#------------------------------------------------------------------------------------------------------------------------------------------

my ($doc1a,$doc2a,$doc3a,$doc4a);
my $row=-1;
my $charge1;
my $charge2;
foreach my $id1 (@list_A){	
	
	my $res1=$doc1->SubUnits->Item($id1)->Name;
	my $chain1=$doc1->SubUnits->Item($id1)->Ancestors->Chain->Name;
	
	foreach my $id2 (@list_U){	
	
		my $res2=$doc1->SubUnits->Item($id2)->Name;
		my $chain2=$doc1->SubUnits->Item($id2)->Ancestors->Chain->Name;
		my($dist,$atom1,$atom2)=res_dist($doc1,$id1,$id2);
	
		if($dist>5 and $dist<=8){				
			$row=$row+1;
			
			print "$res1,$res2,$dist \n";

			$charge1=formal_charge($doc1,$id1);
			$charge2=formal_charge($doc1,$id2);
			
			my @caps1=res_caps($doc1,$id1,$lim1_A,$lim2_A,1);
			my @caps2=res_caps($doc1,$id2,$lim1_U,$lim2_U,1);
	
			my $c1=-1;
			my $c2=-1;
			
			if(substr($res1,0,3) eq "CYS"){$c1=cap_extra($doc1,$id1)}
			if(substr($res2,0,3) eq "CYS"){$c2=cap_extra($doc1,$id2)}
			
			if($c1>=0){
				push(@caps1,$c1);
				print "$id1,$res1,@caps1 \n";
			}
			
			if($c2>=0){	
				push(@caps2,$c2);	
				print "$id2,$res2,@caps2 \n";
			}
			
			($doc1a,$doc2a,$doc3a,$doc4a)=create_mfcc_files_cysteine($doc1,$id1,$id2,$dist,"U1");
			mfcc_cysteine($doc1,$id1,$id2,\@caps1,\@caps2); 
				
			my @list_water=water_list($id1,$id2,\@caps1,\@caps2,\@water_array);
	  		copy_water($doc,\@list_water);
		
			$table2->cell($row,0) = substr($res1,0,index($res1, '_'));
			$table2->cell($row,1) = substr($res2,0,index($res2, '_'));
			$table2->cell($row,2) = $dist;
			$table2->cell($row,3) = $chain1;
			$table2->cell($row,4) = $atom1;
			$table2->cell($row,5) = $atom2;
			$table2->cell($row,6) = $chain2; 
			$table2->cell($row,7) = $charge1;
			$table2->cell($row,8) = $charge2;
	
		}
	}
}


#------------------------------------------------------------------------------------------------------------------------------------------
#                                                   Fim do programa principal/Inicio das funcoes
#------------------------------------------------------------------------------------------------------------------------------------------


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
	
	for ($i=1;$i<=$caps_number;$i++){
		
#            -------------------------
#                   Caps à direita
#            -------------------------
		
		my $id_cap_r = $id_residue + $i;
		
		my $ok="false";
		
		my $cap_right="lim_up";
		if($id_cap_r <= $lim_up){$cap_right=$doc_aux->SubUnits->Item($id_cap_r)->Name}
		
		my $C=$doc_aux->SubUnits->Item($id_residue)->Atoms->DisplayRange("C");
		my $nbonds=$C->AttachedAtoms->Count;
		
		for(my $j=0;$j<=$nbonds-1;$j++){
			my $res2=$C->AttachedAtoms->Item($j)->Ancestors->SubUnit->Name;
			if($res2 eq $cap_right){$ok="true"}
		
		}
		
		
		if($id_cap_r <= $lim_up and $ok eq "true"){push(@caps,$id_cap_r)}
	
#            -------------------------
#                   Caps à esquerda
#            -------------------------
			
		my $id_cap_l = $id_residue - $i;

		my $ok="false";
		
		my $cap_left="lim_down";
		if($id_cap_l >= $lim_down){$cap_left=$doc_aux->SubUnits->Item($id_cap_l)->Name}
		
		my $N=$doc_aux->SubUnits->Item($id_residue)->Atoms->DisplayRange("N");
		my $nbonds=$N->AttachedAtoms->Count;
	
		for(my $j=0;$j<=$nbonds-1;$j++){
			my $res2=$N->AttachedAtoms->Item($j)->Ancestors->SubUnit->Name;
			if($res2 eq $cap_left){$ok="true"}
		}

		if($id_cap_l >= $lim_down and $ok eq "true"){push(@caps,$id_cap_l)}
	
	}
	
	return @caps;
}

#-------------------------------------------------------------------------
#                           sub res_caps
#-------------------------------------------------------------------------

sub cap_extra{
	
	my $doc_aux=@_[0];
	my $id1=@_[1];
	
	my $cap_extra=-1;
	
#            -------------------------
#                     Cap extra
#            -------------------------
		
	my $res= $doc_aux->SubUnits->Item($id1)->Name;	
	my $SG1=$doc_aux->SubUnits->Item($id1)->Atoms->DisplayRange("SG");
	my $nbonds=$SG1->AttachedAtoms->Count;

		for(my $i=0;$i<=$nbonds-1;$i++){
			my $res2=$SG1->AttachedAtoms->Item($i)->Ancestors->SubUnit->Name;
			my $SG2=$SG1->AttachedAtoms->Item($i)->Name;
		
			if($SG2 eq "SG"){
				$cap_extra=id_residue($doc_aux,$res2);
			}	
		}

	return $cap_extra;
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
#                           sub res_dist_files
#-------------------------------------------------------------------------

sub res_dist_files{

	my $doc_aux=@_[0];
	my $doc2_aux=@_[1];
	my $i=@_[2];
	my $j=@_[3];
	
	my $dist_min=1000;
	my $atom1_name;
	my $atom2_name;
	
	foreach my $atom1 (@{$doc_aux->SubUnits->Item($i)->Atoms}){
			foreach my $atom2 (@{$doc2_aux->SubUnits->Item($j)->Atoms}){
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
#                      sub create_mfcc_files_cysteine
#-------------------------------------------------------------------------

sub create_mfcc_files_cysteine{

		my $doc_aux=@_[0];
		my $id1=@_[1];
		my $id2=@_[2];
		my $dist=@_[3];
		my $aux=@_[4];
		
		my $dist=int($dist);
		
		my $res1=$doc_aux->SubUnits->Item($id1)->Name;
		my $res2=$doc_aux->SubUnits->Item($id2)->Name;

		$res1=substr($res1,0,index($res1, '_'));
		$res2=substr($res2,0,index($res2, '_'));
		
		my $name;
		if($dist < 10){$name="00"}                            
		if($dist > 9  and $dist < 100){$name="0"}  
		if($dist > 99 and $dist < 1000){$name=""}
		
		my $new_name=$aux.'_D_'.$name.$dist . '_R_'.$res1.'_'.$res2;
	    my $new_name2= $new_name . "_A";		
		my $doc1a=Documents->New("$new_name2.xsd");	
	
		$new_name2= $new_name . "_B";
		my $doc2a=Documents->New("$new_name2.xsd");
		
		$new_name2= $new_name . "_C";
		my $doc3a=Documents->New("$new_name2.xsd");
		
	    $new_name2= $new_name . "_D";
		my $doc4a=Documents->New("$new_name2.xsd");

		return($doc1a,$doc2a,$doc3a,$doc4a);
}

#-------------------------------------------------------------------------
#                      sub change_res_name
#-------------------------------------------------------------------------

sub change_res_name{

	my $doc_aux=@_[0];
	
	foreach my $residue (@{$doc_aux->SubUnits}){
		if(substr($residue->Name,0,3) ne "HOH"){
			my $chain=$residue->Ancestors->Chain->Name;
			$residue->Name=$residue->Name.'_'.$chain;
		}
	}
	
	return $doc_aux;
}

#-------------------------------------------------------------------------
#                           sub id_residue
#-------------------------------------------------------------------------

sub id_residue{
	
	my $doc_aux=@_[0];
	my $residue_name=@_[1];
	
	my $id_residue;
	my $id_aux=-1;
	foreach my $residue (@{$doc_aux->SubUnits}){
		$id_aux=$id_aux+1;
		my $name = $residue->Name;
		if($residue_name eq "$name"){
			$id_residue=$id_aux;
			last
		}
	}

	return $id_residue;	
}

#-------------------------------------------------------------------------
#                           sub mfcc_cysteine
#-------------------------------------------------------------------------

sub mfcc_cysteine{ 
	
	my $doc_aux=@_[0];
	my $id1=@_[1];
	my $id2=@_[2];
	my @caps1=@{$_[3]};
	my @caps2=@{$_[4]};

	my $res1=$doc_aux->SubUnits->Item($id1)->Name;
	my $res2=$doc_aux->SubUnits->Item($id2)->Name;
	
	my @id_list;
	
	push(@id_list,$id1);
	push(@id_list,$id2);
	push(@id_list,@caps1);
	push(@id_list,@caps2);
	
	#--------------------------------------------------------------------
	#                   Arquivo A  ->  k*-R1-k + C*-R2-C    
	#--------------------------------------------------------------------

	$doc1a->CopyFrom($doc_aux);
	
	my $id=-1;
	foreach my $residue (@{$doc1a->SubUnits}){
		$id=$id+1;
		my $ok="true";
		
		foreach my $id_aux (@id_list){	
			if($id == $id_aux){
				$ok="false";
				last;
			}
		}
		if($ok eq "true"){$doc1a->SubUnits->DisplayRange($residue->Name)->Delete}
	}

	foreach my $atom (@{$doc1a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}
	
	#--------------------------------------------------------------------
	#                  Arquivo B  ->   k*-R1-k +  C*- C
	#--------------------------------------------------------------------
	
	$doc2a->CopyFrom($doc1a);
	$doc2a->SubUnits->DisplayRange($res2)->Delete;

	foreach my $atom (@{$doc2a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}

	#--------------------------------------------------------------------
	#                  Arquivo C  ->  k*-k + C*-R2-C  
	#--------------------------------------------------------------------
	
	$doc3a->CopyFrom($doc1a);
	$doc3a->SubUnits->DisplayRange($res1)->Delete;
	
	foreach my $atom (@{$doc3a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}

	#--------------------------------------------------------------------
	#                       Arquivo D  ->   k*-k + C*-C   
	#--------------------------------------------------------------------
	
	$doc4a->CopyFrom($doc1a);
	$doc4a->SubUnits->DisplayRange($res1)->Delete;
	$doc4a->SubUnits->DisplayRange($res2)->Delete;

	foreach my $atom (@{$doc4a->Atoms}){
		if ($atom->Name eq "N" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "C" and $atom->NumBonds<3){$atom->AdjustHydrogen}
		if ($atom->Name eq "SG" and $atom->NumBonds<2){$atom->AdjustHydrogen}
	}
	
}

#-------------------------------------------------------------------------
#                           sub print_2d
#-------------------------------------------------------------------------

sub print_2d {
	my @array_2d=@_;
	for(my $i = 0; $i <= $#array_2d; $i++){
	   for(my $j = 0; $j <= $#{$array_2d[$i]} ; $j++){
	      print "$array_2d[$i][$j] ";
	   }
	   print "\n";
	}
}

#-------------------------------------------------------------------------
#                          sub length_cube
#-------------------------------------------------------------------------

sub length_cube{
	
	my $doc_aux=@_[0];
		
	my $dmax=-10000;
	foreach my $atom (@{$doc_aux->Atoms}){
		my $x=$atom->X;
		my $y=$atom->Y;
		my $z=$atom->Z;
	
		if($x>$dmax){$dmax=abs($x)}
		if($y>$dmax){$dmax=abs($y)}
		if($z>$dmax){$dmax=abs($z)}
	
	}

	return $dmax;	
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
	
	my $id_cube;
	my $idx;
	my $idy;
	my $idz;
	
	my @cube_array;
	
	my $id_residue=-1;
	foreach my $residue (@{$doc_aux->SubUnits}){
		$id_residue=$id_residue+1;
		
		my $name=$residue->Name;
		
		if(substr($name,0,3) eq "HOH"){
		
			my $atom=$residue->Atoms->Item(0); 
								
			my $x=$atom->X;
			my $y=$atom->Y;
			my $z=$atom->Z;
					
			$idx=int($x/$r)+1;
			$idy=int($y/$r)+1;
			$idz=int($z/$r)+1;
			
			$id_cube = $idx + $L*($idy-1) + $L*$L*($idz-1);
			
			push(@{$cube_array[$id_cube]},$id_residue);
			
			#printf "%10d %10d %10d %10d %10.2f %10.2f %10.2f %10d %10d %10d \n",$r,$L,$id,$id_atom,$x,$y,$z,$idx,$idy,$idz;
			#print "$id_cube,$id_residue,$name \n";
		}
		
	}
	return @cube_array;
}

#-------------------------------------------------------------------------
#                       sub pdb_quadrante1
#-------------------------------------------------------------------------

sub pdb_quadrante1{
	
	my $doc_aux=@_[0];
	
	my $x;
	my $y;
	my $z;
	
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

	$doc_aux->CreateSet("all_atoms", \@atom_set); 
	$doc_aux->Sets("all_atoms")->Translate(Point(X => $dx, Y => $dy, Z => $dz));

	return $doc_aux;
}

#-------------------------------------------------------------------------
#                          sub water_list
#-------------------------------------------------------------------------

sub water_array{
	
	my $doc_aux=@_[0];      #pdb sem água
	my $doc2_aux=@_[1];     #pdb completo
	my $r_cube=@_[2];
	my $L=@_[3];
	my @cube_array=@{$_[4]};
	
	my $N=$L*$L;
	my @water_array;
	
	my $id_residue=-1;
	foreach my $residue (@{$doc_aux->SubUnits}){
		$id_residue=$id_residue+1;
			
		my $center = $residue->CenterOfGeometry;
		my $x=$center->X;
		my $y=$center->Y;
		my $z=$center->Z;
		
		my $id=id_cube($r_cube,$L,$x,$y,$z);
	
		my @v;
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
	
		foreach my $cube (@v){	
			if($cube >=1 and $cube <= ($L*$L*$L)){
				my $nwater=$#{$cube_array[$cube]};
				for(my $i = 0; $i <= $nwater; $i++){
					my $id_water=$cube_array[$cube][$i];
					my($dist,$atom1,$atom2)=res_dist_files($doc_aux,$doc2_aux,$id_residue,$id_water);								
					if($dist<=2.5){push(@{$water_array[$id_residue]},$id_water)}					
				}				
			}
		}
	}

	return @water_array;
}	
	
#--------------------------------------------------------------------
#              			     sub uniq   
#--------------------------------------------------------------------
	
sub uniq {
    my %seen;
    grep !$seen{$_}++, @_;
}

#--------------------------------------------------------------------
#              			   sub water_list   
#--------------------------------------------------------------------

sub water_list{
	  
	my $id1=$_[0];
	my $id2=$_[1];
	my @caps1=@{$_[2]};
	my @caps2=@{$_[3]};
	my @water_array=@{$_[4]};

	my @list_water;
		
	if($#{$water_array[$id1]} >= 0){push(@list_water,@{$water_array[$id1]})}					
	if($#{$water_array[$id2]} >= 0){push(@list_water,@{$water_array[$id2]})}	
	
	foreach my $id_cap1 (@caps1){
		if($#{$water_array[$id_cap1]} >= 0){push(@list_water,@{$water_array[$id_cap1]})};
	}
	
	foreach my $id_cap2 (@caps2){
		if($#{$water_array[$id_cap2]} >= 0){push(@list_water,@{$water_array[$id_cap2]})};
	}
	
	@list_water=uniq(@list_water);
	
	return @list_water;
}

#--------------------------------------------------------------------
#              			   sub copy_water   
#--------------------------------------------------------------------

sub copy_water{
	my $doc_aux=@_[0];
	my @water_list=@{$_[1]};

	foreach my $id_water (@water_list){	
		$doc1a->CopyFrom($doc_aux->SubUnits->Item($id_water));
		$doc2a->CopyFrom($doc_aux->SubUnits->Item($id_water));
		$doc3a->CopyFrom($doc_aux->SubUnits->Item($id_water));
		$doc4a->CopyFrom($doc_aux->SubUnits->Item($id_water));
	}
}

#--------------------------------------------------------------------
#              			   sub formal_charge   
#--------------------------------------------------------------------

sub formal_charge{
	my $doc_aux=@_[0];
	my $id_residue=@_[1];

	my $formalChargeSum=0;
	foreach my $atom (@{$doc_aux->SubUnits->Item($id_residue)->Atoms}){
		my $formalCharge = $atom->FormalCharge->Value;
		$formalChargeSum += $formalCharge;
	}
	
	return $formalChargeSum;
}




