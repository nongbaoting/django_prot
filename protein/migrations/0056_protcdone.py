# Generated by Django 2.2.5 on 2022-04-24 00:06

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0055_scopecla'),
    ]

    operations = [
        migrations.CreateModel(
            name='protCDOne',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('protin_id', models.CharField(db_index=True, max_length=200)),
                ('cdd_nameCat', models.TextField()),
                ('cdd_idCat', models.TextField()),
                ('protin_nr', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='nrinfo_protCDOne', to='protein.NrInfo')),
            ],
        ),
    ]
