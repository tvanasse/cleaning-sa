function combinations = snpm_enumerate_combinations(num_subjects,num_group)
    num_combinations  = nchoosek(num_subjects, num_group);
    combinations      = zeros(num_combinations, num_subjects);
    combination       = 1:num_subjects;
    
    for combination_num = 1:num_combinations
        group_2_index = num_group + 1;
        subject_num   = 1;
        
        for group_1_index = 1:num_group
            while subject_num < combination(group_1_index)
                combination(group_2_index) = subject_num;
                group_2_index              = group_2_index + 1;
                subject_num                = subject_num + 1;
            end
            
            subject_num = combination(group_1_index) + 1;
        end
        
        for group_2_index = group_2_index:num_subjects
            combination(group_2_index) = group_2_index;
        end
        
        combinations(combination_num, :) = combination;
        
        for index = num_group:-1:1
            combination(index) = combination(index) + 1;
            
            if combination(index) <= num_group + index
                for broken_index = (index + 1):num_group
                    combination(broken_index) = combination(broken_index - 1) + 1;
                end
                
                break;
            end
        end
    end
end
