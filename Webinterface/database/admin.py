from django.contrib import admin

# Register your models here.
from .models import Gene, Technology

class TechnologyInline (admin.StackedInline):
    model = Technology
    fields = ['technology', 'motif_file', 'snp_file', 'gene']
    extra = 1
 
class TechnologyAdmin(admin.ModelAdmin):   
    fields = ['technology', 'motif_file', 'gene']
    
class GeneAdmin(admin.ModelAdmin):
    fieldsets = [
                 ('Name', {'fields': ['gene_name', 'refFile']}),
                 ]
    inlines = [TechnologyInline]
    
admin.site.register(Gene, GeneAdmin)
admin.site.register(Technology, TechnologyAdmin)