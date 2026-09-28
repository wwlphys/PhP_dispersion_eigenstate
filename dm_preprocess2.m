fid_in = fopen(qpoints.yaml, 'r');
fid_out = fopen(clean_qpoints.txt, 'w');

if fid_in == -1 || fid_out == -1
    error('文件打开失败');
end

fprintf('preprocessing qpoints.yaml...\n');
line_count = 0;

while ~feof(fid_in)
    line = fgetl(fid_in);
    line_count = line_count + 1;
    
    if mod(line_count, 100000) == 0
        fprintf('进度: 已处理 %d 行\n', line_count);
    end
    
    if ischar(line)
        % 使用正则表达式提取所有数字（包括负数）
        numbers = regexp(line, '-?\d+\.?\d*', 'match');
        
        if ~isempty(numbers)
            % 将提取的数字字符串写入文件
            for i = 1:length(numbers)
                fprintf(fid_out, '%s ', numbers{i});
            end
            fprintf(fid_out, '\n');
        end
    end
end

fclose(fid_in);
fclose(fid_out);
fprintf('处理完成！共处理 %d 行\n', line_count);
