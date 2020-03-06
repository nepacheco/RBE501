function T = get_fwdkin(dh_table,is_sym)
    rows = size(dh_table,1);
    if is_sym
        T = sym('T',[4,4,rows]);
    else
        T = zeros(4,4,rows);
    end
    for i = 1:rows
        if i == 1
            T(:,:,i) = tdh(dh_table(i,:));
        else
            T(:,:,i) = T(:,:,i-1)*tdh(dh_table(i,:));
            if is_sym
                T(:,:,i) = simplify(T(:,:,i),'Steps',20);
            end
        end
    end
end