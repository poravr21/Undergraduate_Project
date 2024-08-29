clc;
close all;
clear all;
%%%SIMULATION PARAMETERS
NOISE_VAR   = 2.511*(10^(-13)) ;
NUM_PILOTS  = 4 ; %20 to 10
AP_ANTENNAS = 1 ;
side     = .250 ;
NUM_APs     = 64;
NUM_USERS   = 12 ;  % 20 to 40
p_Transmit  = .300 ;
CHANNEL_REALIZATIONS = 100;
NUM_Schemes = 1 ;
delta       = 0.7 ; %% for AP selection
MIN_INTER   = 0 ;
Rho_u       = p_Transmit/NOISE_VAR ;    %Normalized signal-to-noise ratio
Bar_4 = zeros(2,NUM_Schemes) ;
SETUPS=5;
initial=0;
for each_setup = 1:SETUPS
    % Intialization
    R_corr          = zeros(AP_ANTENNAS,AP_ANTENNAS,NUM_APs,NUM_USERS) ;
    H               = zeros(NUM_APs,NUM_USERS,CHANNEL_REALIZATIONS) ;
    ntkl            = zeros(NUM_APs,NUM_PILOTS,CHANNEL_REALIZATIONS) ;
    % AP Selection
    [A , beta] = Function_ApSelection(NUM_APs, NUM_USERS, side, delta);
    
    %% SPATIAL CORRELATION MATRIX
    R_corr = Spatial_Correlation_Matrix(A, NUM_USERS,NUM_APs,AP_ANTENNAS);
    random_pilot_assignment = Function_random_pilot(NUM_USERS,NUM_PILOTS ); %% Random pilot assignment Scheme.
    %% Channel Generation
    for k = 1:NUM_USERS
        for l = find(A(:,k))'
            %applying large scale fading coefficient
            H(l,k,:) = sqrt(beta(l,k)/2)*(randn(1,1,CHANNEL_REALIZATIONS)+sqrt(-1)*randn(1,1,CHANNEL_REALIZATIONS));
        end
    end
    %% AWGN channel
    ntkl = sqrt(NOISE_VAR/2)*(randn(NUM_APs,NUM_PILOTS,CHANNEL_REALIZATIONS)+sqrt(-1)*randn(NUM_APs,NUM_PILOTS,CHANNEL_REALIZATIONS));
    
    %% CHANNEL ESTIMATION
    Hhat_mmse = zeros(NUM_APs*AP_ANTENNAS,NUM_USERS,CHANNEL_REALIZATIONS);
    NUM_NMSE = zeros(NUM_APs, NUM_USERS,CHANNEL_REALIZATIONS);
    DEN_NMSE = zeros(NUM_APs, NUM_USERS,CHANNEL_REALIZATIONS);
    NMSE_sample_mmse = zeros(1,CHANNEL_REALIZATIONS);
    Gammmma_each = zeros(NUM_APs*AP_ANTENNAS,NUM_USERS,CHANNEL_REALIZATIONS);
    Gammaa       = zeros(NUM_APs*AP_ANTENNAS,NUM_USERS);
   
   %% MMSE channel estimation 
    for each_sam = 1:CHANNEL_REALIZATIONS
        
        for l = 1:NUM_APs
            
            for t = 1:NUM_PILOTS
                
                ztkl           = sqrt(p_Transmit*NUM_PILOTS)*sum(H(l,t==random_pilot_assignment,each_sam),2)+ntkl(l,t,each_sam); %%eta
                Psi_tkl        = (p_Transmit*NUM_PILOTS*sum(R_corr(:,:,l,t==random_pilot_assignment),4)+eye(AP_ANTENNAS)); %%eta
                
                % channel estimation of users that use pilot t
                for k = find(t==random_pilot_assignment)
                    
                    % only for those existing channels
                    b_linear      = R_corr(:,:,l,k)/Psi_tkl;
                    Hhat_mmse(l,k,each_sam) = sqrt(p_Transmit*NUM_PILOTS)*b_linear*ztkl;  %MMSE esimate Hhat=b_linear*Rx_signal %%eta
                end
            end
        end  
        %% finding Normalized mean square error (NMSE)
        for k = 1:NUM_USERS
            for l = find(A(:,k))' 
                NUM_NMSE(l,k,each_sam)  = abs(H(l,k,each_sam)-Hhat_mmse(l,k,each_sam)).^2;
                DEN_NMSE(l,k,each_sam)  = abs(H(l,k,each_sam)).^2;
            end
        end
        NMSE_sample_mmse(each_sam)      = sum(sum(NUM_NMSE(:,:,each_sam),2),1)./sum(sum(DEN_NMSE(:,:,each_sam),2),1);
    end
    NMSE_mmse_Scheme(each_setup)   = sum(NMSE_sample_mmse)/(CHANNEL_REALIZATIONS);  %% NMSE for each setup
end
Avg_NMSE = sum(NMSE_mmse_Scheme)/SETUPS; %% average NMSE
figure(1);
bar(NMSE_mmse_Scheme)
disp('The average NMSE =')
disp(Avg_NMSE)
