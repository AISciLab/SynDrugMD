import pandas as pd

df1 = pd.read_csv("ALL_Data_withGene1_withID_S.csv")
df2 = pd.read_csv("Human_SL.csv")
df1_expanded = df1.copy()
df1_expanded['Anchor Target'] = df1_expanded['Anchor Target'].apply(lambda x: x.split(';'))
df1_expanded['library Target'] = df1_expanded['library Target'].apply(lambda x: x.split(';'))

# 将df2转换为字典
df2_dict = {(row['n1.identifier'], row['n2.identifier']): True for _, row in df2.iterrows()}
df2_dict.update({(row['n2.identifier'], row['n1.identifier']): True for _, row in df2.iterrows()})


# 检查df1中的成分对是否在df2_dict中
def check_pairs(row, pairs_dict):
    a_all = row['Anchor Target']
    b_all = row['library Target']
    print(row['id'])
    for a_ in a_all:
        for b_ in b_all:
            if a_.isdigit() and b_.isdigit():
                if pairs_dict.get((int(a_), int(b_)), False):
                    return 'y'
    return 'n'


# 添加匹配结果列
df1['match'] = df1_expanded.apply(check_pairs, axis=1, args=(df2_dict,))
df1.to_csv("All_Data_withGene1_withID_S_match.csv", index=False)

# 统计match列中n的行数
n_count = df1['match'].value_counts().get('y', 0)

print(df1)
print(f"满足条件的行数: {n_count}")
