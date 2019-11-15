## This little function compares column names between two given datasets
## (data1 and data2), checking which variables of data1 are NOT in data2

## Just a very lazy function to avoid typing multiple parenthesis, brackets etc

checkols <- function(data1, data2) {
        diffs <- colnames(data1)[which(!colnames(data1) %in% colnames(data2))]
        print(diffs)
}