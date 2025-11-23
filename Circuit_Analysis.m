set(0, 'DefaultFigureVisible', 'off');
%% we are using modified nodal analysis to form equations , then we use matlab symbolic math solver to solve the formed equations(also used in formiing equations), then ode15i to solve the differential equation
%%delete previous simulation result to save new one
if(exist('Results.txt'))
    delete('Results.txt');
end
if(exist('Voltages_graph.fig'))
    delete('Voltages_graph.fig');
end
if(exist('Currents_graph.fig'))
    delete('Currents_graph.fig');
end
if(exist('Cap_curr.fig'))
    delete('Cap_curr.fig');
end
if(exist('Res_C.fig'))
    delete('Res_C.fig');
end

%we take and read input here
%% --- FIX 1: Hardcoded filename to match Flask backend ---
fname='circuit_netlist.txt';
netlist_file=fopen(fname);
%% ---------------------------------------------------------
netlist=textscan(netlist_file,'%s %s %s %s');
fclose(netlist_file);
fileID1=fopen('Element_indep.txt','wt+'); % file from which we will read the values , node , element type
%Initialize
num_Elements=0; %no of elements
num_V=0; %no of voltage sources
num_I=0; %no of current sources
num_Nodes=0; %no of nodes, excluding ground (node 0)
num_L=0; %no of inductors
num_C=0;%no of capacitor
num_R=0;%no of resistor
%%in short the follwing code is parsing the netlist for us to process
for i=1:length(netlist{1})
    s=netlist{1}{i};
    switch(s(1))
        case{'R','L','C','V','I'}
            fprintf(fileID1,[netlist{1}{i} ' ' netlist{2}{i} ' ' netlist{3}{i} ' ' netlist{4}{i} '\n']);
    end
end
%%reading the netlist from Element_indep.txt file
[Name,N1,N2,value]=textread('Element_indep.txt','%s %s %s %s');
for i=1:length(Name)
    switch(Name{i}(1))
        case{'R','L','C'}
            num_Elements=num_Elements+1;
            Element(num_Elements).Name=Name{i};
            Element(num_Elements).Node1=str2num(N1{i});
            Element(num_Elements).Node2=str2num(N2{i});
            Element(num_Elements).Value=str2double(value{i});
            if(Name{i}(1)=='L')
                num_L=num_L+1;
                Inductor(num_L).Name=Name{i};
                Inductor(num_L).N1=str2num(N1{i});
                Inductor(num_L).N2=str2num(N2{i});
                Inductor(num_L).Value=str2double(value{i});
            end
            if(Name{i}(1)=='R')
                num_R=num_R+1;
                Resistor(num_R).Name=Name{i};
                Resistor(num_R).N1=str2num(N1{i});
                Resistor(num_R).N2=str2num(N2{i});
                Resistor(num_R).Value=str2double(value{i});
            end
            if(Name{i}(1)=='C')
                num_C=num_C+1;
                Capacitor(num_C).Name=Name{i};
                Capacitor(num_C).N1=str2num(N1{i});
                Capacitor(num_C).N2=str2num(N2{i});
                Capacitor(num_C).Value=str2double(value{i});
            end
        case{'V'}
            num_V=num_V+1;
            Volt_source(num_V).Name=Name{i};
            Volt_source(num_V).Node1=str2num(N1{i});
            Volt_source(num_V).Node2=str2num(N2{i});
            Volt_source(num_V).Value=str2double(value{i});

        case{'I'}
            num_I=num_I+1;
            Current_source(num_I).Name=Name{i};
            Current_source(num_I).Node1=str2num(N1{i});
            Current_source(num_I).Node2=str2num(N2{i});
            Current_source(num_I).Value=str2double(value{i});
    end
    num_Nodes=max(str2num(N1{i}),max(str2num(N2{i}),num_Nodes));
end

%%close and delete the file
fclose(fileID1);
delete('Element_indep.txt');

%%creating the equations for the independent voltage sources and applying KCL at the nodes

node_equation=cell(num_Nodes,1); 
volt_equation=cell(num_V,1);
for i=1:num_V
    switch((Volt_source(i).Node1==0)||(Volt_source(i).Node2==0))
        case{1}% if one of the node is ground
            if(Volt_source(i).Node1==0)
                volt=['v_' num2str(Volt_source(i).Node2) '=' '-' num2str(Volt_source(i).Value)];
                node_equation{Volt_source(i).Node2}=[node_equation{Volt_source(i).Node2} '-' 'i_' Volt_source(i).Name];
            else
                volt=['v_' num2str(Volt_source(i).Node1) '='  num2str(Volt_source(i).Value)];
                node_equation{Volt_source(i).Node1}=[node_equation{Volt_source(i).Node1} '+' 'i_' Volt_source(i).Name];
            end
            volt_equation{i}=volt;
        case{0}%if none of the node is ground
            volt=['v_' num2str(Volt_source(i).Node1) '-' 'v_' num2str(Volt_source(i).Node2) '=' num2str(Volt_source(i).Value)];
            volt_equation{i}=volt;
            node_equation{Volt_source(i).Node1}=[node_equation{Volt_source(i).Node1} '+' 'i_' Volt_source(i).Name];
            node_equation{Volt_source(i).Node2}=[node_equation{Volt_source(i).Node2} '-' 'i_' Volt_source(i).Name];
    end
end
%A flag used for deciding which solver to finally use , 0  meaning purely resistive circuiit , 1 with rc, rl, rlc ,lc
solver_flag=0;
%add the passive element currents using KCL to the node equations, and make the equations for inductors
L_equation=cell(num_L,1);
L_ctr=0;
for i=1:num_Elements
    switch(Element(i).Name(1))
        case{'R'}
            switch((Element(i).Node1==0)||(Element(i).Node2==0))
                case{0}
                    node_equation{Element(i).Node1}=[node_equation{Element(i).Node1} '+' '(' 'v_' num2str(Element(i).Node2) '-' 'v_' num2str(Element(i).Node1) ')' '/' num2str(Element(i).Value)];
                    node_equation{Element(i).Node2}=[node_equation{Element(i).Node2} '+' '(' 'v_' num2str(Element(i).Node1) '-' 'v_' num2str(Element(i).Node2) ')' '/' num2str(Element(i).Value)];
                case{1}
                    if(Element(i).Node1==0)
                        node_equation{Element(i).Node2}=[node_equation{Element(i).Node2} '-' '(' 'v_' num2str(Element(i).Node2) ')' '/' num2str(Element(i).Value)];
                    else
                        node_equation{Element(i).Node1}=[node_equation{Element(i).Node1} '-' '(' 'v_' num2str(Element(i).Node1) ')' '/' num2str(Element(i).Value)];
                    end
            end
        case{'C'}
            if(solver_flag==0)
                solver_flag=1;
            end
            switch((Element(i).Node1==0)||(Element(i).Node2==0))
                case{0}
                    node_equation{Element(i).Node1}=[node_equation{Element(i).Node1} '+' num2str(Element(i).Value) '*' '(' 'vp(' num2str(Element(i).Node2) ')' '-' 'vp(' num2str(Element(i).Node1) ')' ')'];
                    node_equation{Element(i).Node2}=[node_equation{Element(i).Node2} '+' num2str(Element(i).Value) '*' '(' 'vp(' num2str(Element(i).Node1) ')' '-' 'vp(' num2str(Element(i).Node2) ')' ')'];
                case{1}
                    if(Element(i).Node1==0)
                        node_equation{Element(i).Node2}=[node_equation{Element(i).Node2} '-' num2str(Element(i).Value) '*' '(' 'vp(' num2str(Element(i).Node2) ')' ')'];
                    else
                        node_equation{Element(i).Node1}=[node_equation{Element(i).Node1} '-' num2str(Element(i).Value) '*' '(' 'vp(' num2str(Element(i).Node1) ')' ')'];
                    end
            end
        case{'L'}
            if(solver_flag==0)
                solver_flag=1;
            end
            L_ctr=L_ctr+1;
            switch((Element(i).Node1==0)||(Element(i).Node2==0))
                case{0}
                    node_equation{Element(i).Node1}=[node_equation{Element(i).Node1} '-' 'i_' Element(i).Name];
                    node_equation{Element(i).Node2}=[node_equation{Element(i).Node2} '+' 'i_' Element(i).Name];
                    L_equation{L_ctr}=['v_' num2str(Element(i).Node1) '-' 'v_' num2str(Element(i).Node2) '-' '('  num2str(Element(i).Value) '*' 'ip(' num2str(L_ctr) ')' ')'];
                case{1}
                    if(Element(i).Node1==0)
                        node_equation{Element(i).Node2}=[node_equation{Element(i).Node2}  '+' 'i_' Element(i).Name];
                        L_equation{L_ctr}=['-' 'v_' num2str(Element(i).Node2) '-' '('  num2str(Element(i).Value) '*' 'ip(' num2str(L_ctr) ')' ')'];
                    else
                        node_equation{Element(i).Node1}=[node_equation{Element(i).Node1}  '-' 'i_' Element(i).Name];
                        L_equation{L_ctr}=['v_' num2str(Element(i).Node1) '-' '('  num2str(Element(i).Value) '*' 'ip(' num2str(L_ctr) ')' ')'];
                    end
            end
    end
end
%%Add the independent current sources using KCL to the node equations
for i=1:num_I
    switch((Current_source(i).Node1==0)||(Current_source(i).Node2==0))
        case{1}
            if(Current_source(i).Node1==0)
                node_equation{Current_source(i).Node2}=[node_equation{Current_source(i).Node2} '+' num2str(Current_source(i).Value)];
            else
                node_equation{Current_source(i).Node1}=[node_equation{Current_source(i).Node1} '-' num2str(Current_source(i).Value)];
            end
        case{0}
            node_equation{Current_source(i).Node1}=[node_equation{Current_source(i).Node1} '-' num2str(Current_source(i).Value)];
            node_equation{Current_source(i).Node2}=[node_equation{Current_source(i).Node2} '+' num2str(Current_source(i).Value)];
    end
end
%if solver flag=0 (purely resistive circuit), add the RHS('=0') to each

if(solver_flag==0)
    for i=1:length(node_equation)
        if ~isempty(node_equation{i})
            node_equation{i}=[node_equation{i} '=' '0'];
        end
    end

elseif(solver_flag==1)
    for i=1:num_Nodes %For each nodal KCL equation (only LHS)
        if ~isempty(node_equation{i})
            for j=1:num_Nodes
                node_equation{i}=strrep(node_equation{i},['v_' num2str(j)],['v(' num2str(j) ')']);
            end
            for j=1:num_V
                node_equation{i}=strrep(node_equation{i},['i_' Volt_source(j).Name],['v(' num2str(num_Nodes+j) ')']);
            end
            for j=1:num_L
                node_equation{i}=strrep(node_equation{i},['i_' Inductor(j).Name],['v(' num2str(num_Nodes+num_V+j) ')']);
            end
        end
    end
    for i=1:num_V %For each independent voltage source equation
        for j=1:num_Nodes
            volt_equation{i}=strrep(volt_equation{i},['v_' num2str(j)],['v(' num2str(j) ')']);
        end
        volt_equation{i}=strrep(volt_equation{i},'=','-'); %Modify each independent voltage source equation to only LHS [no RHS ('=0')]
    end


    for i=1:num_L %For each inductor equation (only LHS)
        for j=1:num_Nodes
            L_equation{i}=strrep(L_equation{i},['v_' num2str(j)],['v(' num2str(j) ')']);
        end
    end
end
eqn = cell(num_Nodes + num_V + num_L, 1);

for i = 1:num_Nodes
    if ~isempty(node_equation{i})
        eqn{i} = str2sym(node_equation{i});
    else
        eqn{i} = sym(0);
    end
end

for i = 1:num_V
    eqn{num_Nodes + i} = str2sym(volt_equation{i});
end

for i = 1:num_L
    eqn{num_Nodes + num_V + i} = str2sym(L_equation{i});
end

switch(solver_flag)
    case{0}
        %Create the symbolic variables for node voltages and currents through voltage sources
        variables='syms';
        for i=1:num_Nodes
            variables=[variables ' ' 'v_' num2str(i)];
        end
        for i=1:num_V
            variables=[variables ' ' 'i_' Volt_source(i).Name];
        end
        eval(variables);
        %creating a row vector var of the symbolic variables created above - to be used in solve
        var_string=['var=[' variables(6:end) ']'];
        eval(var_string);
        %Create the symbolic variables for the symbolic equations
        equations='syms';
        for i=1:(num_Nodes+num_V)
            equations=[equations ' ' 'eqn' num2str(i)];
        end
        eval(equations);
        %Create a row vector eqn_solve of the equation symbolic variables
        interm_string=['eqn_solve=[' equations(6:end) ']'];
        eval(interm_string);
        %Assign the equation symbolic variables with the corresponding symbolic equations
        for i=1:(num_Nodes+num_V)
             if ~isempty(eqn{i})
                eqn_string=['eqn' num2str(i) '=' 'eqn{' num2str(i) '}'];
                eval(eqn_string);
             else
                eqn_string=['eqn' num2str(i) '=' '0'];
                eval(eqn_string);
             end
        end
        %Solve the symbolic linear equations using solve
        sol=solve(eval(eqn_solve),var);
        %Note :- We use eval(eqn_solve) to substitute the symbolic equation associated with
        %each equation symbolic variable
        F=fopen('Results.txt','wt+'); %Create an empty text file Results.txt
        fprintf(F,['File name : ' fname]);
        fprintf(F,'\n');
        fprintf(F,'DC ANALYSIS RESULTS \n');
        fprintf(F,'-------------------\n');
        fprintf(F,'NODE VOLTAGES \n');
        for i=1:num_Nodes
            fprintf(F,['v_' num2str(i) ' = ']);
            fprintf(F,num2str(eval(eval(['sol.v_' num2str(i)]))));
            fprintf(F,' V\n');
        end
        if(num_V~=0)
            fprintf(F,'CURRENTS THROUGH INDEPENDENT VOLTAGE SOURCES (NEGATIVE TO POSITIVE TERMINAL) \n');
            for i=1:num_V
                fprintf(F,['i_' Volt_source(i).Name ' = ']);
                fprintf(F,num2str(eval(eval(['sol.i_' Volt_source(i).Name]))));
                fprintf(F,' A\n');
            end
        end
        %  currents through each resistor
        fprintf(F, 'CURRENTS THROUGH RESISTORS \n');
        for i = 1:num_R
            % Assuming Resistor(i).N1 and Resistor(i).N2 represent the nodes connected to the resistor
            if(Resistor(i).N1~=0)
                V_N1 = eval(['sol.v_' num2str(Resistor(i).N1)]); % Voltage at node N1
            else
                V_N1 = 0;
            end
            if(Resistor(i).N2~=0)
                V_N2 = eval(['sol.v_' num2str(Resistor(i).N2)]); % Voltage at node N2
            else
                V_N2 = 0;
            end

            % Current through the resistor using Ohm's Law: I = (V_N1 - V_N2) / R
                current_through_resistor = (V_N1 - V_N2) / Resistor(i).Value;

                fprintf(F, ['I through ' Resistor(i).Name ' = ']);
                fprintf(F, num2str(double(current_through_resistor)));
                fprintf(F, ' A\n');
            end

            %% --- FIX 4: Removed redundant type('Results.txt') ---
            fclose(F); %Close the Results.txt text file

        case{1}
            %Create the state variables for node voltages, currents through voltage sources and inductor currents
            variables='syms';
            for i=1:(num_Nodes+num_V+num_L)
                variables=[variables ' ' 'v' num2str(i) '(t)'];
            end
            eval(variables);
            %Create a row vector var of the state variables - to be used in daeFunction
            var_string=['var=[' variables(6:end) ']'];
            eval(var_string);
            %Convert the symbolic equations (only LHS) to a form suitable for daeFunction
            %Use the converted symbolic equations to make a row vector eqn_daeFunction - to be used in daeFunction
            eqn_string='eqn_daeFunction=[';
            for i=1:length(eqn)
                interm_string=char(eqn{i});
                for j=1:(num_Nodes+num_V+num_L)
                    interm_string=strrep(interm_string,['v(' num2str(j) ')'],['v' num2str(j) '(t)']);
                end
                for j=1:num_Nodes
                    interm_string=strrep(interm_string,['vp(' num2str(j) ')'],['diff(v' num2str(j) '(t)' ',t)']);
                end
                for j=1:num_L
                    interm_string=strrep(interm_string,['ip(' num2str(j) ')'],['diff(v' num2str(num_Nodes+num_V+j) '(t)' ',t)']);
                end
                eqn_string=[eqn_string interm_string ','];
            end
            eqn_string=[eqn_string ']'];
            eval(eqn_string);
            %Use daeFunction to create the function handle odefun
            odefun=daeFunction(eqn_daeFunction,var);

            %% --- THIS IS THE NEW SECTION TO FIX THE SOLVER ---
            v0_guess=zeros(length(eqn_daeFunction),1); % GUESS for initial conditions v
            vp0_guess=zeros(length(eqn_daeFunction),1); % GUESS for initial conditions v'

        % Give a better guess for the voltage source to help decic
        for i=1:num_V
            if(Volt_source(i).Node1 ~= 0 && Volt_source(i).Node2 == 0)
               v0_guess(Volt_source(i).Node1) = Volt_source(i).Value;
            elseif(Volt_source(i).Node2 ~= 0 && Volt_source(i).Node1 == 0)
               v0_guess(Volt_source(i).Node2) = -Volt_source(i).Value;
            end
        end

        %% --- FIX 2: Removed the input prompt for tf. ---
        %% The variable 'tf' is now passed from the Flask backend.
        disp(['Transient analysis from t=0 to t=' num2str(tf)]);
        %% ------------------------------------------------

        options=odeset('RelTol',1e-03,'AbsTol',1e-03);

        % --- NEW FIX: Use decic to find consistent initial conditions ---
        [v0, vp0] = decic(odefun, 0, v0_guess, [], vp0_guess, []);
        % ---------------------------------------------------------

        [t,v]=ode15i(odefun,[0 tf],v0,vp0,options);
        %% --- END OF NEW SECTION ---


        % Calculate currents through resistors and capacitors
        resistor_currents = zeros(length(t), num_R);
        if(num_R~=0)
            for i = 1:num_R

                % Calculate current for resistor: I = (V_node1 - V_node2) / R
                if(Resistor(i).N1~=0)
                    nodeR1_voltage = v(:, Resistor(i).N1:Resistor(i).N1);
                else
                    nodeR1_voltage = zeros(length(t),1);
                end

                if(Resistor(i).N2~=0)
                    nodeR2_voltage = v(:, Resistor(i).N2:Resistor(i).N2);
                else
                    nodeR2_voltage = zeros(length(t),1);
                end
                resistor_currents(:, i) = (nodeR1_voltage - nodeR2_voltage) / Resistor(i).Value;
            end
            figure;
            plot(t,resistor_currents(:,1:num_R));
            Labels = cell(1, num_R);
            for k = 1:num_R
                Labels{k} = Resistor(k).Name;
            end

            legend(Labels);
            title('Current through resistors');
            xlabel('TIME (s)');
            ylabel('CURRENTS (A)');
            exportgraphics(gcf, 'Res_C.png', 'Resolution', 150, 'BackgroundColor', 'black');
            end

        capacitor_currents = zeros(length(t),num_C);
        if(num_C~=0)
            for i = 1:num_C
                if(Capacitor(i).N1~=0)
                    node1_voltage = v(:, Capacitor(i).N1);
                else
                    node1_voltage = zeros(length(t),1);
                end

                if(Capacitor(i).N2~=0)
                    node2_voltage = v(:, Capacitor(i).N2);
                else
                    node2_voltage = zeros(length(t),1);
                end

                voltage_diff = node1_voltage - node2_voltage;

                % Calculate derivative using finite differences (except for first point)
                dvdt = zeros(length(t),1);
                dvdt(2:end) = diff(voltage_diff) ./ diff(t);
                dvdt(1) = dvdt(2); % Approximate first point

                capacitor_currents(:,i) = Capacitor(i).Value * dvdt;
            end
            figure;
            plot(t,capacitor_currents(:,1:num_C));
            legendLabels = cell(1, num_C);
            for k = 1:num_C
                legendLabels{k} = Capacitor(k).Name;
            end

            legend(legendLabels);
            title('Current through Capacitors');
            xlabel('TIME (s)');
            ylabel('CURRENTS (A)');
            exportgraphics(gcf, 'Cap_curr.png', 'Resolution', 150, 'BackgroundColor', 'black');
        end

        table_heading=cell(1,(1+num_Nodes+num_V+num_L));
        table_heading{1}='Time';
        for j=1:num_Nodes
            table_heading{1+j}=['v_' num2str(j)];
        end
        for j=1:num_V
            table_heading{1+num_Nodes+j}=['i_' Volt_source(j).Name];
        end

        for j=1:num_L
            table_heading{1+num_Nodes+num_V+j}=['i_' Inductor(j).Name];
        end
        figure;
        plot(t,v(:,1:num_Nodes)); %Plot the node voltages vs. time
        legend_voltage='legend(';
        for i=1:num_Nodes
            interm_string=table_heading{1+i};
            interm_string=strrep(interm_string,'_','\_');
            legend_voltage=[legend_voltage '''' interm_string '''' ','];
        end
        legend_voltage(end)=')';
        eval(legend_voltage);
        title('Node voltages')
        xlabel('TIME (s)');
        ylabel('NODE VOLTAGES (V)');
        exportgraphics(gcf, 'Voltages_graph.png', 'Resolution', 150, 'BackgroundColor', 'black');

        if((num_V~=0)||(num_L~=0))
            figure; %Create new figure window
            plot(t,v(:,(num_Nodes+1):end)); %Plot the currents through voltage sources and inductor currents vs. time
            legend_current='legend(';
            for i=1:num_V
                interm_string=table_heading{1+num_Nodes+i};
                interm_string=strrep(interm_string,'_','\_');
                legend_current=[legend_current '''' interm_string '''' ','];
            end

            for i=1:num_L
                interm_string=table_heading{1+num_Nodes+num_V+i};
                interm_string=strrep(interm_string,'_','\_');
                legend_current=[legend_current '''' interm_string '''' ','];
            end
            legend_current(end)=')';
            eval(legend_current);
            title('Current through inductor and voltage source')
            xlabel('TIME (s)');
            ylabel('CURRENTS (A)');
            exportgraphics(gcf, 'Currents_graph.png', 'Resolution', 150, 'BackgroundColor', 'black');
        end

        %% --- FIX 3: Write final transient values to Results.txt ---
        F=fopen('Results.txt','wt+');
        fprintf(F,['File name : ' fname]);
        fprintf(F,'\n');
        fprintf(F,'TRANSIENT ANALYSIS - FINAL VALUES at t = %.4f s\n', tf);
        fprintf(F,'-------------------------------------------------\n');

        fprintf(F,'NODE VOLTAGES \n');
        for i=1:num_Nodes
            fprintf(F, '%s = %.5f V\n', table_heading{i+1}, v(end,i));
        end

        if(num_V~=0)
            fprintf(F,'CURRENTS THROUGH INDEPENDENT VOLTAGE SOURCES (NEGATIVE TO POSITIVE TERMINAL) \n');
            for i=1:num_V
                 fprintf(F, '%s = %.5f A\n', table_heading{1 + num_Nodes + i}, v(end, num_Nodes + i));
            end
        end

        if(num_L~=0)
            fprintf(F,'CURRENTS THROUGH INDUCTORS \n');
             for i=1:num_L
                 fprintf(F, '%s = %.5f A\n', table_heading{1 + num_Nodes + num_V + i}, v(end, num_Nodes + num_V + i));
            end
        end

        if(num_R~=0)
            fprintf(F, 'CURRENTS THROUGH RESISTORS \n');
            for i = 1:num_R
                fprintf(F, 'I through %s = %.5f A\n', Resistor(i).Name, resistor_currents(end, i));
            end
        end

        if(num_C~=0)
            fprintf(F, 'CURRENTS THROUGH CAPACITORS \n');
            for i = 1:num_C
                fprintf(F, 'I through %s = %.5f A\n', Capacitor(i).Name, capacitor_currents(end, i));
            end
        end

        fclose(F);
        %% -------------------------------------------------------------
end