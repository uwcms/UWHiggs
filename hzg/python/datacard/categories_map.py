
#this maps human readable category types to
#names in the ttree and a RooFit categories string

categories_map = {'det_based_4cat':
                  {'leafname':'r94cat',
                   'catstring':'r94cat[cat1=1,cat2=2,cat3=3,cat4=4]',
                   'categories':[1,2,3,4],
                   'bkg_models':{'RooBernstein':{1:4, #category:order
                                                 2:5,
                                                 3:5,
                                                 4:5}
                                 }
                   }
                  }

#returns the skeleton of the factory string
def make_background_for_cat(cattype,cat,bkg_model):
    available_cats = categories_map[cattype]['categories']
    available_models = categories_map[cattype]['bkg_models']
    if cat in available_cats:
        if bkg_model in available_models:
            model_info = available_models[bkg_model]
            factory_string = '%s::background_model_cat_%i'%(bkg_model,cat)
            factory_string += '(%s)'
            specialize = ''
            if bkg_model == 'RooBernstein':
                nconsts = model_info[cat]
                pconsts = ['bkg_cat_%i_p%i'%(cat,i+1) for i in range(nconsts)]
                specialize += 'Mzg,1.0,%s'%(','.join(pconst))
            #fill our factory string with useful information
            factory_string = factory_string%specialize
            return factory_string                
        else:
            raise Exception('model %s not in %s!'%(bkg_model,cattype))
    else:
        raise Exception('category %s not in %s!'%(cat,cattype))
