ConfusionMatrix := proc(predicted::{list,rtable,DataSeries}, 
                        actual::{list,rtable,DataSeries}, 
                        {mapping::list:=NULL}, 
                        {difference::truefalse:=false}) #Undocumented option
    local CM, CM_results, CM_split, i, j, labels, n, unique;

    if numelems(predicted) <> numelems(actual) then
        error "predicted must be of same length as actual";
    elif numelems(convert(predicted,'set')) > numelems(convert(actual,'set')) then
        error "unidentified value in prediction";
    end if:

    CM := Array(1..numelems(predicted), 1..4); 
    #Original Values
    CM[..,1] := convert(actual, 'Array');
    #Cluster Predictions
    CM[..,2] := convert(predicted, 'Array'); 

    unique := [op(convert(CM[..,1], 'set'))];

    if mapping <> NULL then
        n := numelems(unique);
        if hastype(mapping, 'list'({'algebraic','string'})) then
            #Original Values mapped onto Cluster Numbers
            CM[..,3] := subs(seq(mapping[i] = i, i = 1..n), CM[..,1]);
            labels := mapping;
        elif hastype(mapping, 'list'(`=`)) then
            #Original Values mapped onto Cluster Numbers
            CM[..,3] := subs(seq(i, i in mapping), CM[..,1]);
            labels := map(lhs, mapping);
        end if;
    else
        n := numelems(unique);
        CM[..,3] := CM[..,1];
        labels := map(convert,unique,'string');
    end if:
    unique := [op(convert(CM[..,3], 'set'))];

    if difference then
        for i from 1 to numelems(predicted) do
            if CM[i,3] = CM[i,2] then
                CM[i,4] := true;
            else
                CM[i,4] := false;
            end if;
        end do;
    end if;

    CM_split := Statistics:-SplitByColumn(CM, 3);

    CM_results := Matrix(n, n, 0):
    for i to n do
        for j to n do
            CM_results[i,j] := numelems(select['flatten'](x -> x = unique[j], CM_split[i][.., 2]));  
        end do:
    end do:
    
    if difference = false then
        return DataFrame(CM_results, 'columns'=labels, 'rows'=labels);
    else
        return CM[..,4];
    end if;
end proc:

ConfusionTable := proc(confusionmatrix::DataFrame, 
                       variable::{string,name},
                      {summary::truefalse:=false})
    local accuracy, CT, dv, DOR, FDR, FNR, FOR, FPR, LRN, LRP, NPV, PPV, prevalence, TNR, TPR, total;

    if ColumnLabels(confusionmatrix) <> RowLabels(confusionmatrix) then
        error "DataFrame is not symmetric";
    end if;

    dv := proc(x,y)
        try
            x/y;
        catch:
            return 0;
        end try;
    end proc:

    CT := Matrix(1..2,1..2):
    #True Positive
    CT[1,1] := confusionmatrix[variable,variable];
    #False Negatives
    CT[1,2] := add(confusionmatrix[variable,..]) - CT[1,1];
    #False Positives
    CT[2,1] := add(confusionmatrix[..,variable]) - CT[1,1];
    #True Negatives
    CT[2,2] := add(add(confusionmatrix[i,j], i=1..upperbound(confusionmatrix)[1]), j=1..upperbound(confusionmatrix)[2]) - CT[1,1] - CT[1,2] - CT[2,1];
    total := CT[1,1] + CT[1,2] + CT[2,1] + CT[2,2];
    accuracy := (CT[1,1] + CT[2,2])/total;
    prevalence := (CT[1,2] + CT[2,1])/total;
    #Sensitivity - True Positive Rate
    TPR := dv(CT[1,1],(CT[1,1] + CT[1,2]));
    #False Negative Rate
    FNR := dv(CT[1,2],(CT[1,1] + CT[1,2]));
    #False Positive Rate
    FPR := dv(CT[2,1],(CT[2,1] + CT[2,2]));
    #Specificity - True Negative Rate
    TNR := dv(CT[2,2],(CT[2,1] + CT[2,2]));
    #Precision - Positive Predictive Value
    PPV := dv(CT[1,1],(CT[1,1] + CT[2,1]));
    #False Discoverty Rate
    FDR := dv(CT[2,1],(CT[1,1] + CT[2,1]));
    #False Omission Rate
    FOR := dv(CT[1,2],(CT[1,2] + CT[2,2]));
    #Negative Predictive Value
    NPV := dv(CT[2,2],(CT[2,1] + CT[2,2]));
    #Positive Likelihood Ratio
    LRP := dv(TPR,FPR);
    #Negative Likelihood Ratio
    LRN := dv(FNR,TNR);
    #Diagnostic Odds Ratio
    DOR := dv(LRP,LRN);
    DocumentTools:-Tabulate(Matrix([[sprintf("%d true positives", CT[1,1]), 
                                     sprintf("%d false negatives", CT[1,2])], 
                                    [sprintf("%d false positives", CT[2,1]), 
                                     sprintf("%d true negatives", CT[2,2])]]), 
                            'width'=50, 'widthmode'='percentage');
    if summary then   
        printf("Accuracy: %.2f\n", accuracy);
        printf("Prevalence: %.2f\n", prevalence);

        printf("Sensitivity: %.2f\n", TPR);
        printf("False Negative Rate: %.2f\n", FNR);
        printf("False Positive Rate: %.2f\n", FPR);
        printf("Specificity: %.2f\n", TNR);

        printf("Positive Predictive Value: %.2f\n", PPV);
        printf("False Discovery Rate: %.2f\n", FDR);
        printf("False Omission Rate: %.2f\n", FOR);
        printf("Negative Predictive Value: %.2f\n", NPV);

        printf("Positive Likelihood Ratio: %.2f\n", LRP);
        printf("Negative Likelihood Ratio: %.2f\n", LRN);
        printf("Diagnostic Odds Ratio: %.2f", DOR);
    end if;
    return NULL;
end proc:
