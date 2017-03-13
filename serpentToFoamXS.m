clear all
clc


fprintf(' \n serpentToFoamXS - Version 5/3/2013\n\nThis script converts a SERPENT(2) output to the crossSections dictionary for neutronDiffusionFoam\n\nWARNING: use Matlab!\n\n');

mev2j = 1.602176487e-13;

cm2m = 0.01;

%caricare il file
rad = input('\nSerpent case name\ne.g.: msfr_run_1\n', 's');
resfilestr = strcat(rad,'_res');
if exist(resfilestr)
        fprintf('\nreading %s.m\n',resfilestr);
        eval(resfilestr)
else
        fprintf('\nThe file %s.m does not exists.\n\nSTOP\n',resfilestr);
        break
end

coreState = input('\n Core state you want the xs to be prepared for:\n N for nominal, \n R for radially expanded, \n A for axially expanded, \n T for Doppler broadened, \n C for expanded coolant, \n CL for expanded cladding  \n', 's');


%creare il file crossSections
if (strcmp('N',coreState))
	globalfilestr = strcat('nuclearData');
elseif (strcmp('R',coreState))
	globalfilestr = strcat('nuclearDataRadialExp');
elseif (strcmp('A',coreState))
	globalfilestr = strcat('nuclearDataAxialExp');
elseif (strcmp('T',coreState))
	globalfilestr = strcat('nuclearDataFuelTemp');
elseif (strcmp('C',coreState))
	globalfilestr = strcat('nuclearDataRhoCool');
elseif (strcmp('CL',coreState))
	globalfilestr = strcat('nuclearDataCladExp');
end

if exist(globalfilestr)
        yesno = input('\nWarning: il file verra sovrascritto\ncontinuo? (yes/no)\n', 's');
        if ((strcmp('yes',yesno)==0))
                fprintf('\nSTOP (no file modified).\n');
                break
        end
end


%beta_eff energy correction
fprintf('\nTreatment of delayed neutron:\nDo you want to use the physical fractions and the real spectrum or the effective fractions and the prompt neutron spectrum?\n');
zeroeff = input('\n zero for physical fractions, eff for energy corrected fractions\n', 's');

if (((strcmp('zero',zeroeff)==0) || (strcmp('eff',zeroeff))))
        zeroeff = input('\nPlease choose:\n zero for physical delayed neutrons fraction and energy spectrum\n eff for energy corrected fractions and CHI_delayed = CHI_prompt\n', 's');
end

if (((strcmp('zero',zeroeff)==0) || (strcmp('eff',zeroeff))))
        fprintf('\n Sorry: invalid %s option.\n',zeroeff);
        break
end

fprintf('\nOpening file %s\n',globalfilestr);
fid=fopen(globalfilestr,'w');

fprintf(fid,'\n/*\n crossSection dictionary\n Generatd by serpentToFoamXS \n %s\n From SERPENT results file: %s\n*/\n',date,resfilestr);

if (strcmp('zero',zeroeff))
        fprintf(fid,'\n/*\nphysical delayed neutron fraction and spectrum\n*/\n');
end

if (strcmp('eff',zeroeff))
        fprintf(fid,'\n/*\neffective delayed neutron fraction, prompt neutron spectrum\n*/\n');
end


fprintf(fid,'\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    location    constant;\n    object      %s;\n}\n',resfilestr);

%numero gruppi
ng = GC_NE(idx);
fprintf('\nNumber of energy groups: %i\n', ng);

%numero ritardati
%nd = PRECURSOR_GROUPS(idx);
nd = length (FWD_ANA_BETA_ZERO)/2 -1;
fprintf('\nNumber of delayed neutron precursor group: %i\n', nd);



fprintf(fid,'\n energyGroups %i ; \n',ng);
fprintf(fid,'\n precGroups %i ; \n',nd);

if (strcmp('N',coreState))
	pTarget = input('\n Target power for eigenvalue calculations:\n');
	keff = input('\n keff:\n');
	fprintf(fid,'\n pTarget %.6e ; \n',pTarget);
	fprintf(fid,'\n keff %.6e ; \n \n',keff);
end
if (strcmp('R',coreState))
	expansionFromNominalR = input('\n Relative radial expansion compared to nominal:\n');
	fprintf(fid,'\n expansionFromNominal %.6e ; \n',expansionFromNominalR);
	radialOrientationX = input('\n Orientation of radial direction, x component:\n');
	radialOrientationY = input('\n Orientation of radial direction, y component:\n');
	radialOrientationZ = input('\n Orientation of radial direction, z component:\n');
	fprintf(fid,'\n radialOrientation %i %i %i ; \n',radialOrientationX,radialOrientationY,radialOrientationZ);
	AxialOrientationX = input('\n Orientation of Axial direction, x component:\n');
	AxialOrientationY = input('\n Orientation of Axial direction, y component:\n');
	AxialOrientationZ = input('\n Orientation of Axial direction, z component:\n');
	fprintf(fid,'\n axialOrientation %i %i %i ; \n',AxialOrientationX,AxialOrientationY,AxialOrientationZ);

end
if (strcmp('A',coreState))
	expansionFromNominalA = input('\n Relative axial expansion compared to nominal:\n');
	fprintf(fid,'\n expansionFromNominal %.6e ; \n',expansionFromNominalA);
end

% Data specific different perturbed states
if (strcmp('T',coreState))
	Tfuelref = input('\n Tfuelref:\n');
	fprintf(fid,' Tfuelref %.6e ; \n',Tfuelref);
	TfuelPerturbed = input('\n TfuelPerturbed:\n');
	fprintf(fid,' TfuelPerturbed %.6e ; \n',TfuelPerturbed);
end
if (strcmp('C',coreState))
	rhoCoolRef = input('rhoCoolRef:\n');
	fprintf(fid,' rhoCoolRef %.6e ; \n',rhoCoolRef);
	rhoCoolPerturbed = input('\n rhoCoolPerturbed:\n');
	fprintf(fid,' rhoCoolPerturbed %.6e ; \n',rhoCoolPerturbed);
end
if (strcmp('CL',coreState))
	Tcladref = input('Tcladref:\n');
	fprintf(fid,' Tcladref %.6e ; \n',Tcladref);
	TcladPerturbed = input('\n TcladPerturbed:\n');
	fprintf(fid,' TcladPerturbed %.6e ; \n',TcladPerturbed);
end
%loop over materials
fprintf(fid,'zones \n ( \n');


NumberOfUniverses = size(GC_UNIVERSE_NAME,1)
stop = false;
while(stop == false)
	universeNumber = input('\n Number of the universe you want to include\n e.g.: 0\n STOP to stop \n', 's');


	if(strcmp(universeNumber,'STOP'))
		stop = true;
		break;
	end

	universeFound = false;
	for i = 1:NumberOfUniverses
		if(strcmp(strcat(GC_UNIVERSE_NAME(i,:)),universeNumber))
			idx=i;
			universeFound = true;
		end
	end
	if (universeFound == false)
		fprintf('\n Warning: universe not found, moving to the next one \n');
	else
		universeName = input('\n Name of this universe in OpenFOAM\n e.g.: fuelIn\n', 's');
		fprintf(fid,'%s \n{ \n',universeName);

		if (strcmp('N',coreState))
			fuelFraction = input('\n volumetric fraction fuel(necessary to evaluate vol power IN the fuel)\n e.g.:  0.3 \n');
			fprintf(fid,' fuelFraction %.6e ; \n',fuelFraction);
		end

		% 1 su V
		fprintf(fid,' IV nonuniform List<scalar> %i (',ng);

		IV = zeros(ng,1);
		for i = 1:ng
				IV(i) = RECIPVEL(idx,(i*2)+1);
				fprintf(fid,'%.6e ',IV(i)/cm2m);
		end
		fprintf(fid,' );');

		%coeff diffusione
		fprintf(fid,'\n D nonuniform List<scalar> %i (',ng);

		D = zeros(ng,1);
		for i = 1:ng
				D(i) = INF_DIFFCOEF(idx,(i*2)-1);
				fprintf(fid,'%.6e ',D(i)*cm2m);
		end
		fprintf(fid,' );');

		%nuSigmaEff
		fprintf(fid,'\n nuSigmaEff nonuniform List<scalar> %i (',ng);

		A = zeros(ng,1);
		for i = 1:ng
				A(i) = NSF(idx,(i*2)+1);
				fprintf(fid,'%.6e ',A(i)/cm2m);
		end
		fprintf(fid,' );');

		%sigmaPow
		fprintf(fid,'\n sigmaPow nonuniform List<scalar> %i (',ng);

		A = zeros(ng,1);
		for i = 1:ng
				A(i) = FISSXS(idx,(i*2)+1) * FISSE(idx,(i*2)+1);
				fprintf(fid,'%.6e ',A(i)/cm2m*mev2j);
		end
		fprintf(fid,' );');

		%scattering matrix
		fprintf(fid,'\n scatteringMatrix  %i  %i ( \n',ng,ng);
		MS = zeros(ng,ng);
		for i = 1:ng
				fprintf(fid,' (');
				for j = 1:ng
						MS(i,j) = GPRODXS(idx,2*(j-1)*ng+2*i-1);
						fprintf(fid,' %.6e ',MS(i,j)/cm2m);
				end
				fprintf(fid,')\n ');
		end
		fprintf(fid,');');

		% sigma disappearence (abs+ capture + group transfer below)
		fprintf(fid,'\n sigmaDisapp nonuniform List<scalar> %i (',ng);

		DISAPP = zeros(ng,1);
		for i = 1:ng
				DISAPP(i) = TOTXS(idx,(i*2)+1) - MS(i,i);
				fprintf(fid,'%.6e ',DISAPP(i)/cm2m);
		end
		fprintf(fid,' );');

		%spettro neutroni pronti
		fprintf(fid,'\n chiPrompt nonuniform List<scalar> %i (',ng);
		XP = zeros(ng,1);
		for i = 1:ng
				XP(i) = CHIP(idx,(i*2)-1);
				fprintf(fid,'%.6e ',XP(i));
		end
		fprintf(fid,' );');

		%spettro neutroni ritardati
		fprintf(fid,'\n chiDelayed nonuniform List<scalar> %i (',ng);
		XD = zeros(ng,1);
		for i = 1:ng
				XD(i) = CHID(idx,(i*2)-1);
				if (strcmp('eff',zeroeff))
						fprintf(fid,'%.6e ',XP(i));
				else
						fprintf(fid,'%.6e ',XD(i));
				end
		end
		fprintf(fid,' );');

		%beta
		fprintf(fid,'\n Beta nonuniform List<scalar> %i (',nd);
		BZ = zeros(nd,1);
		BE = zeros(nd,1);
		for i = 1:nd
				BZ(i) = FWD_ANA_BETA_ZERO(idx,(i*2)+1);
				BE(i) = ADJ_IFP_IMP_BETA_EFF(idx,(i*2)+1);
				if (strcmp('eff',zeroeff))
						fprintf(fid,'%.6e ',BE(i));
				else
						fprintf(fid,'%.6e ',BZ(i));
				end
		end
		fprintf(fid,' );');

		%lambda
		fprintf(fid,'\n lambda nonuniform List<scalar> %i (',nd);
		LAM = zeros(nd,1);
		for i = 1:nd
				LAM(i) = FWD_ANA_LAMBDA(idx,(i*2)+1);
				fprintf(fid,'%.6e ',LAM(i));
		end
		fprintf(fid,' );');
		fprintf(fid,'\n } \n');

	end
end

fprintf(fid,'\n ); \n');
fprintf('\nSaving file %s\n',globalfilestr);
fclose(fid);
