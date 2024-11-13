function [output_data, output_data2, output_data3] = ...
    file_data(filename, input_data, input_data2, input_data3)

if nargin==1
    method = 'read_file';
elseif nargin==2 && ischar(input_data)
    method = 'read_folder';
elseif nargin==2
    method = 'single_input_write_file';
else
    method = 'many_structs_write_file';
end
    

switch method

    case 'read_file' % file to data
        data_pool = load(filename);
        if isstruct(data_pool) % single file to single or multiple structs
            switch numel(fieldnames(data_pool))
                case 1
                    output_data = data_pool.data;
                case 2
                    output_data = data_pool.data;
                    output_data2 = data_pool.data2;
                case 3
                    output_data = data_pool.data;
                    output_data2 = data_pool.data2;
                    output_data3 = data_pool.data3;
            end
        else % the content is not a struct, simply output it
            output_data = data_pool;
        end
        
    case 'single_input_write_file' % data to file, no output
    
        if isstruct(input_data) % struct data to file

            data = input_data;
            dt=whos('data');
            if dt.bytes<2e9 % smaller than 2 GB
                save(filename, 'data');
            else
                save(filename, 'data', '-v7.3');
            end

        elseif iscell(input_data) % cell data to file

            fid = fopen(file_name, 'w');
            N_cell = numel(input_data);
            digits = ceil(log10(2*N_cell));
            format = ['%0', num2str(digits), 'd'];

            for i = 1:N_cell
                line = input_data{i};
                for k=1:numel(line)
                     if k>1
                         fprintf(fid, ',');
                     end
                     fprintf(fid, format, line(k));
                 end
                 fprintf(fid, '\n');
            end
            fclose(fid);

        else % matrix data to file
			
            writematrix(input_data, filename);
			

        end
        
    case 'many_structs_write_file'
        
        switch nargin
            
            case 3
                data = input_data;
                data2 = input_data2;
				dt=whos('data');
				dt2=whos('data2');
				if dt.bytes+dt2.bytes<2e9
					save(filename, 'data', 'data2');
				else
					save(filename, 'data', 'data2', '-v7.3');
				end
                
            case 4
                data = input_data;
                data2 = input_data2;
                data3 = input_data3;
				dt=whos('data');
				dt2=whos('data2');
				dt3=whos('data3');
				if dt.bytes+dt2.bytes+dt3.bytes<2e9
					save(filename, 'data', 'data2', 'data3');
				else
					save(filename, 'data', 'data2', 'data3', '-v7.3');
				end
        end
        

    case 'read_folder' % multiple files to a struct
        file_in = dir(filename);
        N_file = length(file_in);
        folder_str = file_in(1).folder;
        
        output_data = struct();

        for i_file = 1:N_file
            name = file_in(i_file).name;
            output_data(i_file).filename = name;
            output_data(i_file).(input_data) = file_data([folder_str, '\', name]);
        end
        
        
    
        
end
    