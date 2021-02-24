function [acc_mat, time_mat] = get_mat_data(data, alp)
    d = data;
    acc_mat = [eval_acc_in_data(d,30,2,alp) eval_acc_in_data(d,30,3,alp) eval_acc_in_data(d,30,4,alp);
               eval_acc_in_data(d,9,2,alp)  eval_acc_in_data(d,9,3,alp) eval_acc_in_data(d,9,4,alp); 
               eval_acc_in_data(d,4,2,alp) eval_acc_in_data(d,4,3,alp) eval_acc_in_data(d,4,4,alp)];
    time_mat = [eval_time_in_data(d,30,2,alp) eval_time_in_data(d,30,3,alp) eval_time_in_data(d,30,4,alp);
               eval_time_in_data(d,9,2,alp)  eval_time_in_data(d,9,3,alp) eval_time_in_data(d,9,4,alp); 
               eval_time_in_data(d,4,2,alp) eval_time_in_data(d,4,3,alp) eval_time_in_data(d,4,4,alp)];
end

function v = eval_acc_in_data(data, D, n, alp)
    d = data;
    acc_str = ['d.accuracy_D' num2str(D) '_n' num2str(n) '_alp' num2str(alp)];
    eval(['v = ' acc_str ';'])
end

function v = eval_time_in_data(data, D, n, alp)
    d = data;
    time_str = ['d.time_D' num2str(D) '_n' num2str(n) '_alp' num2str(alp)];
    eval(['v = ' time_str ';'])
end