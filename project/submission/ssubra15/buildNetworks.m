%Build networks from data using 60 seconds before and after each seizure
%for each event and based on absolute value of pairwise correlation

%List of patients
patient_labels = {'PY04N007','PY04N008','PY04N012','PY04N013','PY04N015','PY05N004','PY05N005','PY11N003','PY11N004','PY11N006'};
num_patients = length(patient_labels);

%Patient outcomes
outcomes = [0;1;0;1;0;0;1;0;0;1]; %1 is a success

%Networks
networks = cell(2,1);

%Label vector initialize
event_labels = [];

seizure = 0; %To count the number of seizures

for p = 1:num_patients
%Load the data
struct_name = ['fsv_pwr' patient_labels{p} '.mat'];
load(struct_name);
struct_name = ['info_' patient_labels{p} '.mat'];
load(struct_name);

    for n = 1:nevents %For each event for that patient
        seizure = seizure + 1; %Increment the number of seizures
        current_event = eval(['snap' num2str(n) '_gamma']); %Pull the current event
        
        %Get number of electrodes
        num_elec = size(current_event,1);

        network = zeros(num_elec,num_elec);

        current_event = current_event(:,start_marks(n)-60:end_marks(n)+60); %Only seizure + 60 seconds before and after
        network = abs(corr(current_event')); %Compute pairwise correlation and abs_val it

        networks{1,seizure} = network;
        networks{2,seizure} = RR_electrodes;
        event_labels = [event_labels; outcomes(p)]; %Update the event label based on patient outcome
    end
end