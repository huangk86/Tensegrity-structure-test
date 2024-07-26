function [EstimatedHarmonics] = ExtractHarmonics (rawSpectrum,HarmonicExtractParameters,DisplayOnOFF)
% Extract the peaks within prescribed frequency ranges 
% The so called phase should be in rad, the so-called angle should be in deg.

fn=HarmonicExtractParameters.fn;
fn_range=HarmonicExtractParameters.fn_range;
Extracted_Harmonics=HarmonicExtractParameters.M; % The number of extracted harmonics is the same with generated harmonics

    Freq=rawSpectrum.Freq;
    Y=abs(rawSpectrum.Y);
    AngY=angle(rawSpectrum.Y)*180/pi;
    F=abs(rawSpectrum.F);
    AngF=angle(rawSpectrum.F)*180/pi;
    
    % find the primary harmonic term from F
    [F_amp,F_index] = max(F );
    F_ang =AngF (F_index);
    
    
    % find harmonics from Y
    for n=1:Extracted_Harmonics
        RangeN=find ((Freq >=(n*fn -fn_range))&(Freq <=(n*fn+fn_range)));
        
        [Y_amp(n),Y_index]=max(Y (RangeN));
        Y_ang(n)=AngY (RangeN(Y_index));       

        
        %         H_amp (n)=Y_amp (n)/(F_amp );
         FreqHarmonics (n)=Freq (RangeN(Y_index));
        %         H_pha (n) =Y_pha (n)-F_pha ;
    end
    
    EstimatedHarmonics.amp=ArrowColumn2Column(Y_amp).';
    EstimatedHarmonics.ang =    (ArrowColumn2Column(Y_ang)-F_ang *ones(length (Y_ang),1)).';  %% The so-called phase shouldnot be a relative phase, -F_pha *ones(length (Y_pha ),1 )
    EstimatedHarmonics.f =       ArrowColumn2Column(F_amp).';
          EstimatedHarmonics.freqharmonics=FreqHarmonics;
     if DisplayOnOFF==1
        disp(['Monitored force amplitude ', num2str(EstimatedHarmonics.f )]);
        disp(['Monitored force phase ', num2str(F_ang ), ' (deg)']);
        disp('  ')
       end   
    %% Validation of the extracted harmonics by visualisation of the results
%     主要看看幅值相位有没有找对。
  if DisplayOnOFF==1
    figure;
    semilogy(Freq ,Y ) ;hold on;
    semilogy(FreqHarmonics ,Y_amp ,'o');hold on;
    xlim([0 2100])
    
    figure
    plot(Freq ,AngY ) ;hold on;
    plot(FreqHarmonics ,Y_ang ,'o');hold on;
    xlim([0 2100])
  end
