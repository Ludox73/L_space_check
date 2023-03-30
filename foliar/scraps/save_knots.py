import taskdb2
import pandas as pd

cols = ['name', 'foliar_tri', 'foliar_orient', 'floer_simple']
dk14 = taskdb2.ExampleDatabase('KnotsInS3').dataframe()[cols]
dk14 = dk14[dk14.name.apply(lambda x: 'n' in x)]
dk15 = taskdb2.ExampleDatabase('KnotsInS3NonAlt15').dataframe()[cols]
dk15.name = dk15.name.apply(lambda x: 'K' + x)
dk16 = taskdb2.ExampleDatabase('KnotsInS3NonAlt16').dataframe()[cols]
dk16.name = dk16.name.apply(lambda x: 'K' + x)

dk = pd.concat([dk14, dk15, dk16]).reset_index()
dk.loc[dk.foliar_tri.isnull(), 'floer_simple'] = 1 
dk.loc[dk.foliar_tri.notnull(), 'floer_simple'] = -1

nonsimp = dk.index[dk.floer_simple==-1]
dk['tets'] = None
dk.loc[nonsimp, 'tets'] = dk.foliar_orient[nonsimp].apply(eval).apply(len)
del dk['id']

dk.to_csv('persistent_knots.csv', index=False)
