%% Defining Constants 
clear;
close all;
bands = [1,2,3,4,8,16,32,64];
filterOrder = [4,4,3,3,3,3,3,3];
freql = 90;
freqh = 5760;
%% Loading/Generating audio files and equilizing their lengths
[sound1, Fs] = audioread('sound1.wav');
[sound2, Fs1] = audioread('sound2.wav');
[sound1,sound2]=equalize_length(sound1,sound2);
[melody1,Fs2]=audioread('melody1.wav');
[melody2,Fs3]=audioread('melody2.wav');
[melody1,melody2]=equalize_length(melody1,melody2);
[sound1d,sound2d] = delay_sounds(sound1,sound2,Fs);
noise = 2*rand(size(sound1,1), 1)-1;
noise=[noise,noise];
%% setting the variables to loop on
sounds1={sound1,noise,sound1,sound1d,melody1};%for envelope
sounds2={noise,sound1,sound2,sound2d,melody2};%for fine structure
sampling_freq=[Fs,Fs,Fs,Fs,Fs2];
chimaers=["speech-noise","noise-speech","speech-speech",...
    "dichotic-chimaera","melody-melody"];
%% main for loop
for j=1:size(chimaers,2)
    if(j==5)
        filterOrder=[4,3,2,2,2,2,2,2];
    end
    sound1=sounds1{1,j};
    sound2=sounds2{1,j};
    f1=figure('name',chimaers(j));
    subplot(4, 1, 1);
    NDFT = 2^nextpow2(size(sound1,1));
    [h,w] = freqz(sound1(:,1), NDFT);
    plot(w/pi,abs(h));
    xlabel('Normalized Freq(pi rad/sample)');
    ylabel('X(w)');
    title('SOUND1 SIGNAL DFT')

    subplot(4, 1, 2);
    t = 0:size(sound1,1)-1;
    plot(t,sound1(:,1));
    xlabel('Time in secs');
    ylabel('Amplitude'); 
    title('SOUND1 SIGNAL TIME DOMAIN')

    subplot(4, 1, 3);
    [h,w] = freqz(sound2(:,2), NDFT);
    plot(w/pi,abs(h));
    xlabel('Normalized Freq(pi rad/sample)');
    ylabel('X(w)');
    title('SOUND2 SIGNAL DFT')

    subplot(4, 1, 4);
    t = 0:size(sound2,1)-1;
    plot(t,sound2(:,2));
    xlabel('Time in secs');
    ylabel('Amplitude'); 
    title('SOUND2 SIGNAL TIME DOMAIN');
    % Please put the path of directory where you want to save the plots
    saveas(f1,"D:\6th sem\EXP5\plots\"+chimaers(j)+"_input.jpg")
    for idx=1:size(bands,2)
        n_bands=bands(idx);
        Fs=sampling_freq(j);
        result=chimaera_synthesizer(sound1,sound2,n_bands,...
        filterOrder(idx),Fs,freqh,freql);
        f2=figure('name',chimaers(j)+"_bands_"+bands(idx));
        NDFT = 2^nextpow2(size(result,1));
        subplot(2,1,1);
        [h,w] = freqz(result(:,1), NDFT);
        plot(w/pi,abs(h));
        xlabel('Normalized Freq(pi rad/sample)');
        ylabel('Y(w)');
        title('PROCESSED AUDIO SIGNA    L DFT')

        subplot(2, 1, 2);
        plot(t,result(:,1));
        xlabel('Time in secs');
        ylabel('Amplitude'); 
        title('PROCESSED AUDIO SIGNAL TIME DOMAIN');

        sgtitle("Result For No of Bands :"+n_bands);
        %please put the directory where you want to save the audio files
        str="D:\6th sem\EXP5\audio_results\";
        str=str++chimaers(j)+"_result_"+n_bands+".wav";
        audiowrite(str, result, Fs);
        %please put the directory where you want to save the plots
        str="D:\6th sem\EXP5\plots\";
        str=str++chimaers(j)+"_result_"+n_bands+".jpg";
        saveas(f2,str);
    end
end
%% function to make chimaeric sound  
function result = chimaera_synthesizer(sound1,sound2,n_bands...
    ,filterOrder,Fs,freqh,freql)
    result=0;
    step=nthroot(freqh/freql,n_bands);
    fl=freql;
    for i= 1:n_bands
        fr=fl*step;
        [b,a] = butter(filterOrder,[2*fl/Fs, 2*fr/Fs], 'bandpass');
        y1 = filter(b, a, sound1); 
        y2= filter(b, a, sound2);
        %SOUND1 using as envelope
        y_hilbert1=hilbert(y1);
        envlp=abs(y_hilbert1);
        %SOUND2 using as fine structure
        y_hilbert2=hilbert(y2);
        fine=y2./abs(y_hilbert2);
        result = result + fine.*envlp;
        fl=fr;
    end
    result = result./max(abs(result));
end
%% function to introduce delay
%Delaying sounds by 700us in opposite directions
function [s1,s2] = delay_sounds(sound1,sound2,Fs)
    delay_time=700E-6; % 700us delay
    dly = ceil(delay_time * Fs);
    
    DataDly2 = zeros(size(sound2,1)+dly, 2);
    q1 = dly:size(sound2,1)+dly-1;
    q2 = 1:size(sound2,1);
    DataDly2(q1,1) = sound2(:,1);%delaying channel 1 of sound2
    DataDly2(q2,2) = sound2(:,2);
    s2=DataDly2; %sound2 delayed
    %Writing delayed sound2 in the sound2d.wav file
    audiowrite("sound2d.wav",s2,Fs); 
  

    DataDly1 = zeros(size(sound1,1)+dly, 2);
    q1 = dly:size(sound1,1)+dly-1;
    q2 = 1:size(sound1,1);
    DataDly1(q1,2) = sound1(:,1);%delaying channel 2 of sound1
    DataDly1(q2,1) = sound1(:,2);
    s1=DataDly1; %sound1 delayed in opposite channel
    audiowrite("sound1d.wav",s1,Fs);
end
%% Function to make the lengths of two sound waves equal
function [s1,s2] = equalize_length(sound1,sound2)
    if(size(sound1,1)~=size(sound2,1))
        if(size(sound1,1)<size(sound2,1))
            sound1=[sound1 ;zeros(size(sound2,1)-size(sound1,1),2)];
        else
            sound2=[sound2;zeros(size(sound1,1)-size(sound2,1),2)];
        end
    end
    s1=sound1;
    s2=sound2;
end